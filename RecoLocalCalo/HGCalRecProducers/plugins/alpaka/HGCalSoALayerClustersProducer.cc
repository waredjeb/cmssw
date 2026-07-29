#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
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

#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalSoALayerClustersProducer : public stream::SynchronizingEDProducer<> {
  public:
    HGCalSoALayerClustersProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          getTokenDeviceRecHits_{consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))},
          getTokenDeviceClusters_{consumes(config.getParameter<edm::InputTag>("hgcalRecHitsLayerClustersSoA"))},
          getTokenMaxLayerPerSide_{consumes<unsigned int>(config.getParameter<edm::InputTag>("hgcalMaxLayerPerSide"))},
          deviceTokenSoAClusters_{produces()},
          thresholdW0_(config.getParameter<double>("thresholdW0")),
          positionDeltaRho2_(config.getParameter<double>("positionDeltaRho2")) {}

    ~HGCalSoALayerClustersProducer() override = default;

    void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
      // Get LayerClusters almost-SoA on device: this has still the same
      // cardinality as the RecHitsSoA, but has all the required information
      // to assemble the clusters, i.e., it has the cluster index assigned to
      // each rechit.
      auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
      auto const inputClusters_v = deviceInputClusters.view();
      //
      // Allocate output SoA for the clusters, one entry for each cluster
      auto device_numclusters = cms::alpakatools::make_device_view<const unsigned int>(
          alpaka::getDev(iEvent.queue()), inputClusters_v.numberOfClustersScalar());
      auto host_numclusters = cms::alpakatools::make_host_view<unsigned int>(num_clusters_);
      alpaka::memcpy(iEvent.queue(), host_numclusters, device_numclusters);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      // Get RecHitsSoA on the device
      auto const& deviceInputRecHits = iEvent.get(getTokenDeviceRecHits_);
      auto const inputRechits_v = deviceInputRecHits.view();

      // Get LayerClusters almost-SoA on device: this has still the same
      // cardinality as the RecHitsSoA, but has all the required information
      // to assemble the clusters, i.e., it has the cluster index assigned to
      // each rechit.
      auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
      auto const inputClusters_v = deviceInputClusters.view();

      // Number of layers per side, as computed by the producer that filled
      // hgcalRecHitsSoA: needed to recover the per-side, z-side-independent
      // `layer` convention (rhtools.getLayerWithOffset(seed)) from the
      // combined (layerOnSide + zside*maxLayerPerSide) value stored per rechit.
      const unsigned int maxLayerPerSide = iEvent.get(getTokenMaxLayerPerSide_);

      reco::CaloClusterDeviceCollection output(
          iEvent.queue(), num_clusters_, num_clusters_, num_clusters_, num_clusters_);
      // Zero-initialise the whole buffer: run() only writes a subset of the
      // columns (position, energy, seedID, algoID, corrected energies), so
      // this guarantees the remaining columns (flags, timing) hold a defined
      // value instead of uninitialised device memory.
      output.zeroInitialise(iEvent.queue());
      auto output_v = output.view();

      // caloID is deliberately left zero-initialised here: it has a virtual
      // destructor, so it is neither constructible on device nor safely
      // copyable to it. HGCalLayerClustersFromSoAProducer sets it host-side
      // when it builds the host cluster collection.

      // Allocate workspace SoA cluster
      HGCalSoAClustersExtraDeviceCollection outputWorkspace(iEvent.queue(), num_clusters_);
      auto output_workspace_v = outputWorkspace.view();

      algo_.run(iEvent.queue(),
                num_clusters_,
                thresholdW0_,
                positionDeltaRho2_,
                maxLayerPerSide,
                inputRechits_v,
                inputClusters_v,
                output_v,
                output_workspace_v);
      iEvent.emplace(deviceTokenSoAClusters_, std::move(output));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("hgcalRecHitsLayerClustersSoA", edm::InputTag("TO BE DEFINED"));
      desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("TO BE DEFINED"));
      desc.add<edm::InputTag>("hgcalMaxLayerPerSide", edm::InputTag("TO BE DEFINED", "maxLayerPerSide"));
      desc.add<double>("thresholdW0", 2.9);
      desc.add<double>("positionDeltaRho2", 1.69);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<HGCalSoARecHitsDeviceCollection> const getTokenDeviceRecHits_;
    device::EDGetToken<HGCalSoARecHitsExtraDeviceCollection> const getTokenDeviceClusters_;
    edm::EDGetTokenT<unsigned int> const getTokenMaxLayerPerSide_;
    device::EDPutToken<reco::CaloClusterDeviceCollection> const deviceTokenSoAClusters_;
    HGCalLayerClustersSoAAlgoWrapper algo_;
    unsigned int num_clusters_;
    float thresholdW0_;
    float positionDeltaRho2_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSoALayerClustersProducer);
