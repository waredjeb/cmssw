#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalSoAClustersExtra.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalTilesConstants.h"

#include "HGCalLayerClustersSoAAlgoWrapper.h"

#include <cstdint>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace {
    // Build a View addressing [offset, offset + size) of a larger cluster
    // collection, so that each subdetector can be filled in place by its own
    // kernel launch instead of being assembled separately and copied.
    ::reco::CaloClusterSoA::View sliceOf(::reco::CaloClusterSoA::View view,
                                         cms::soa::size_type offset,
                                         cms::soa::size_type size) {
      auto position = view.position();
      auto energy = view.energy();
      auto indexes = view.indexes();
      auto timing = view.timing();
      return ::reco::CaloClusterSoA::View{decltype(position){size,
                                                           {position.x().data() + offset},
                                                           {position.y().data() + offset},
                                                           {position.z().data() + offset},
                                                           {position.layer().data() + offset},
                                                           {position.cells().data() + offset}},
                                        decltype(energy){size,
                                                         {energy.energy().data() + offset},
                                                         {energy.correctedEnergy().data() + offset},
                                                         {energy.correctedEnergyUncertainty().data() + offset}},
                                        decltype(indexes){size,
                                                          {indexes.caloID().data() + offset},
                                                          {indexes.algoID().data() + offset},
                                                          {indexes.seedID().data() + offset},
                                                          {indexes.flags().data() + offset}},
                                        decltype(timing){
                                            size, {timing.time().data() + offset}, {timing.timeError().data() + offset}}};
    }
  }  // namespace

  class HGCalSoALayerClustersProducer : public stream::SynchronizingEDProducer<> {
  public:
    HGCalSoALayerClustersProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          deviceTokenSoAClusters_{produces()},
          clusterOffsetsToken_{produces("clusterOffsets")} {
      // One entry per subdetector: the clustering is per-subdetector, but all
      // the clusters end up in a single collection, in the configured order.
      for (auto const& pset : config.getParameter<std::vector<edm::ParameterSet>>("layerClusters")) {
        auto const detector = pset.getParameter<std::string>("detector");
        Input input{consumes(pset.getParameter<edm::InputTag>("hgcalRecHitsSoA")),
                    consumes(pset.getParameter<edm::InputTag>("hgcalRecHitsLayerClustersSoA")),
                    consumes<unsigned int>(pset.getParameter<edm::InputTag>("hgcalMaxLayerPerSide")),
                    static_cast<float>(pset.getParameter<double>("thresholdW0")),
                    static_cast<float>(pset.getParameter<double>("positionDeltaRho2")),
                    detector == "BH",
                    detector == "EE"   ? ::reco::CaloCluster::hgcal_em
                    : detector == "BH" ? ::reco::CaloCluster::hgcal_scintillator
                                       : ::reco::CaloCluster::hgcal_had};
        inputs_.push_back(input);
      }
      numberOfClusters_.resize(inputs_.size(), 0);
    }

    ~HGCalSoALayerClustersProducer() override = default;

    void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
      // The number of clusters of each subdetector is only known on the device.
      // Copy them all back on the same queue, so that the merged output can be sized
      for (size_t d = 0; d < inputs_.size(); ++d) {
        auto const& deviceInputClusters = iEvent.get(inputs_[d].clusters);
        auto device_numclusters = cms::alpakatools::make_device_view<const unsigned int>(
            alpaka::getDev(iEvent.queue()), deviceInputClusters.view().numberOfClustersScalar());
        auto host_numclusters = cms::alpakatools::make_host_view<unsigned int>(numberOfClusters_[d]);
        alpaka::memcpy(iEvent.queue(), host_numclusters, device_numclusters);
      }
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      const size_t numberOfInputs = inputs_.size();

      // Offset of each subdetector in the merged collection; the last entry is
      // the total number of clusters. Published, because the consumers that
      // read the per-subdetector rechit collections need it to map a
      // per-subdetector cluster index into the merged one.
      std::vector<uint32_t> clusterOffsets(numberOfInputs + 1, 0);
      for (size_t d = 0; d < numberOfInputs; ++d) {
        clusterOffsets[d + 1] = clusterOffsets[d] + numberOfClusters_[d];
      }
      const cms::soa::size_type totalNumberOfClusters = clusterOffsets.back();

      reco::CaloClusterDeviceCollection output(
          iEvent.queue(), totalNumberOfClusters, totalNumberOfClusters, totalNumberOfClusters, totalNumberOfClusters);
      // Zero-initialise the whole buffer: run() only writes a subset of the
      // columns (position, energy, seedID, algoID, corrected energies), and it
      // skips a subdetector altogether when that one has no clusters, so this
      // guarantees the remaining rows and columns hold a defined value instead
      // of uninitialised device memory.
      output.zeroInitialise(iEvent.queue());
      auto output_v = output.view();

      // Workspaces are addressed per-subdetector, but cannot be ended
      // until after full completion as asynchronously ran
      std::vector<HGCalSoAClustersExtraDeviceCollection> workspaces;
      workspaces.reserve(numberOfInputs);

      for (size_t d = 0; d < numberOfInputs; ++d) {
        auto const& deviceInputRecHits = iEvent.get(inputs_[d].recHits);
        auto const& deviceInputClusters = iEvent.get(inputs_[d].clusters);

        const unsigned int maxLayerPerSide = iEvent.get(inputs_[d].maxLayerPerSide);

        const cms::soa::size_type size = numberOfClusters_[d];
        workspaces.emplace_back(iEvent.queue(), size);

        algo_.run(iEvent.queue(),
                  size,
                  inputs_[d].thresholdW0,
                  inputs_[d].positionDeltaRho2,
                  maxLayerPerSide,
                  inputs_[d].isScintillator,
                  inputs_[d].algoId,
                  deviceInputRecHits.view(),
                  deviceInputClusters.view(),
                  sliceOf(output_v, clusterOffsets[d], size),
                  workspaces.back().view());
      }

      iEvent.emplace(deviceTokenSoAClusters_, std::move(output));
      iEvent.emplace(clusterOffsetsToken_, std::move(clusterOffsets));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      edm::ParameterSetDescription layerClustersDesc;
      layerClustersDesc.add<edm::InputTag>("hgcalRecHitsLayerClustersSoA", edm::InputTag("TO BE DEFINED"));
      layerClustersDesc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("TO BE DEFINED"));
      layerClustersDesc.add<edm::InputTag>("hgcalMaxLayerPerSide", edm::InputTag("TO BE DEFINED", "maxLayerPerSide"));
      layerClustersDesc.add<double>("thresholdW0", 2.9);
      layerClustersDesc.add<double>("positionDeltaRho2", 1.69);
      layerClustersDesc.ifValue(edm::ParameterDescription<std::string>(
                                    "detector", "EE", true, edm::Comment("the HGCAL component used to create clusters.")),
                                edm::allowedValues<std::string>("EE", "FH", "BH"));
      desc.addVPSet("layerClusters", layerClustersDesc, {});
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    struct Input {
      device::EDGetToken<HGCalSoARecHitsDeviceCollection> recHits;
      device::EDGetToken<HGCalSoARecHitsExtraDeviceCollection> clusters;
      edm::EDGetTokenT<unsigned int> maxLayerPerSide;
      float thresholdW0;
      float positionDeltaRho2;
      bool isScintillator;
      ::reco::CaloCluster::AlgoId algoId;
    };

    std::vector<Input> inputs_;
    device::EDPutToken<reco::CaloClusterDeviceCollection> const deviceTokenSoAClusters_;
    edm::EDPutTokenT<std::vector<uint32_t>> const clusterOffsetsToken_;
    HGCalLayerClustersSoAAlgoWrapper algo_;
    std::vector<unsigned int> numberOfClusters_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSoALayerClustersProducer);
