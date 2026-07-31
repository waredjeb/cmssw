// Rebuild the layer clusters of the HGCAL subdetectors from the portable SoA
// representation.

#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsExtraHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"

#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include <cstdint>
#include <vector>

class HGCalLayerClustersFromSoAProducer : public edm::stream::EDProducer<> {
public:
  HGCalLayerClustersFromSoAProducer(edm::ParameterSet const& config)
      : getTokenSoAClusters_(consumes(config.getParameter<edm::InputTag>("src"))),
        getTokenClusterOffsets_(consumes(config.getParameter<edm::InputTag>("clusterOffsets"))),
        hitsTime_(config.getParameter<unsigned int>("nHitsTime")) {
    for (auto const& pset : config.getParameter<std::vector<edm::ParameterSet>>("layerClusters")) {
      inputs_.push_back({consumes<HGCalSoARecHitsHostCollection>(pset.getParameter<edm::InputTag>("hgcalRecHitsSoA")),
                         consumes<HGCalSoARecHitsExtraHostCollection>(
                             pset.getParameter<edm::InputTag>("hgcalRecHitsLayerClustersSoA")),
                         // timing is not implemented in the digitisation of the scintillator
                         pset.getParameter<std::string>("detector") != "BH"});
    }

    produces<std::vector<float>>("InitialLayerClustersMask");

    produces<reco::CaloClusterHostCollection>();
    produces<ticl::HitsAndFractionsHost>();
  }

  ~HGCalLayerClustersFromSoAProducer() override = default;

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
    const size_t numberOfInputs = inputs_.size();

    auto const& deviceData = iEvent.get(getTokenSoAClusters_);
    auto const deviceView = deviceData.view();
    auto const position_v = deviceView.position();
    auto const energy_v = deviceView.energy();
    auto const indexes_v = deviceView.indexes();
    auto const timing_v = deviceView.timing();
    const int numberOfClusters = position_v.metadata().size();

    // clusterOffsets[d] is the index, in the merged collection, of the first
    // cluster of the d-th subdetector; the last entry is the total number of
    // clusters. It is filled by HGCalSoALayerClustersProducer from the same
    // list of subdetectors, in the same order.
    auto const& clusterOffsets = iEvent.get(getTokenClusterOffsets_);
    if (clusterOffsets.size() != numberOfInputs + 1) {
      throw cms::Exception("Configuration")
          << "HGCalLayerClustersFromSoAProducer is configured with " << numberOfInputs
          << " subdetectors, but the input clusterOffsets describes " << (clusterOffsets.size() - 1)
          << ": the 'layerClusters' list must match the one of the module producing 'src'.";
    }

    // Number of rechits in each cluster, excluding outliers
    std::vector<int> hitsPerCluster(numberOfClusters, 0);
    int numberOfClusteredHits = 0;
    for (size_t d = 0; d < numberOfInputs; ++d) {
      auto const& soaRecHitsExtra_v = iEvent.get(inputs_[d].recHitsExtra).view();
      for (int i = 0; i < soaRecHitsExtra_v.metadata().size(); ++i) {
        const auto clusterIndex = soaRecHitsExtra_v[i].clusterIndex();
        if (clusterIndex == -1) {
          continue;
        }
        ++hitsPerCluster[clusterIndex + clusterOffsets[d]];
        ++numberOfClusteredHits;
      }
    }

    auto clusters = std::make_unique<reco::CaloClusterHostCollection>(
        cms::alpakatools::host(), numberOfClusters, numberOfClusters, numberOfClusters, numberOfClusters);
    auto hitsAndFractions = std::make_unique<ticl::HitsAndFractionsHost>(
        cms::alpakatools::host(), numberOfClusteredHits, numberOfClusters);
    auto clusters_v = clusters->view();
    auto hits_v = hitsAndFractions->view();

    // Prefix-sum the per-cluster counts into the CSR offsets, and copy the
    // per-cluster features across. caloID is the only field that the kernels do
    // not fill; algoID is taken from the SoA, where it is already set per
    // subdetector.
    int hitsOffset = 0;
    for (int i = 0; i < numberOfClusters; ++i) {
      hits_v.offsets()[i].keys_offsets() = hitsOffset;
      hitsOffset += hitsPerCluster[i];

      clusters_v.position()[i].x() = position_v[i].x();
      clusters_v.position()[i].y() = position_v[i].y();
      clusters_v.position()[i].z() = position_v[i].z();
      clusters_v.position()[i].layer() = position_v[i].layer();
      clusters_v.position()[i].cells() = position_v[i].cells();
      clusters_v.energy()[i].energy() = energy_v[i].energy();
      clusters_v.energy()[i].correctedEnergy() = energy_v[i].correctedEnergy();
      clusters_v.energy()[i].correctedEnergyUncertainty() = energy_v[i].correctedEnergyUncertainty();
      clusters_v.indexes()[i].caloID() = reco::CaloID(reco::CaloID::DET_HGCAL_ENDCAP);
      clusters_v.indexes()[i].algoID() = indexes_v[i].algoID();
      clusters_v.indexes()[i].seedID() = indexes_v[i].seedID();
      clusters_v.indexes()[i].flags() = indexes_v[i].flags();
    }
    hits_v.offsets()[numberOfClusters].keys_offsets() = hitsOffset;
    assert(hitsOffset == numberOfClusteredHits);

    // scatter the hits, and collect the per-hit times needed below.
    std::vector<int> cursor(numberOfClusters, 0);
    std::vector<std::vector<float>> times(numberOfClusters);
    std::vector<std::vector<float>> timeErrors(numberOfClusters);
    for (size_t d = 0; d < numberOfInputs; ++d) {
      auto const& soaCells_v = iEvent.get(inputs_[d].cells).view();
      auto const& soaRecHitsExtra_v = iEvent.get(inputs_[d].recHitsExtra).view();
      for (int i = 0; i < soaRecHitsExtra_v.metadata().size(); ++i) {
        const auto clusterIndex = soaRecHitsExtra_v[i].clusterIndex();
        if (clusterIndex == -1) {
          continue;
        }
        const int j = clusterIndex + clusterOffsets[d];
        const auto slot = hits_v.offsets()[j].keys_offsets() + cursor[j]++;
        hits_v.content().values()[slot] = ticl::HitAndFraction{soaCells_v[i].detid(), 1.f};
        if (soaCells_v[i].timeError() < 0.f) {
          continue;
        }
        times[j].push_back(soaCells_v[i].time());
        timeErrors[j].push_back(1.f / (soaCells_v[i].timeError() * soaCells_v[i].timeError()));
      }
    }

    // Assign time to each cluster
    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    for (size_t d = 0; d < numberOfInputs; ++d) {
      for (int j = clusterOffsets[d]; j < static_cast<int>(clusterOffsets[d + 1]); ++j) {
        const auto timeCl = inputs_[d].hasTiming
                                ? timeEstimator.fixSizeHighestDensity(times[j], timeErrors[j], hitsTime_)
                                : std::pair<float, float>(-99.f, -1.f);
        clusters_v.timing()[j].time() = timeCl.first;
        clusters_v.timing()[j].timeError() = timeCl.second;
      }
    }

    auto layerClustersMask = std::make_unique<std::vector<float>>(numberOfClusters, 1.f);

    iEvent.put(std::move(layerClustersMask), "InitialLayerClustersMask");
    iEvent.put(std::move(clusters));
    iEvent.put(std::move(hitsAndFractions));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hltHgcalSoALayerClustersProducer"));
    desc.add<edm::InputTag>("clusterOffsets", edm::InputTag("hltHgcalSoALayerClustersProducer", "clusterOffsets"));
    edm::ParameterSetDescription layerClustersDesc;
    layerClustersDesc.add<edm::InputTag>("hgcalRecHitsLayerClustersSoA",
                                         edm::InputTag("hltHgcalSoARecHitsLayerClustersProducer"));
    layerClustersDesc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("hltHgcalSoARecHitsProducer"));
    layerClustersDesc.ifValue(edm::ParameterDescription<std::string>(
                                  "detector", "EE", true, edm::Comment("the HGCAL component used to create clusters.")),
                              edm::allowedValues<std::string>("EE", "FH", "BH"));
    desc.addVPSet("layerClusters", layerClustersDesc, {});
    desc.add<unsigned int>("nHitsTime", 3);
    desc.add<std::string>("timeClname", "timeLayerCluster");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  struct Input {
    edm::EDGetTokenT<HGCalSoARecHitsHostCollection> cells;
    edm::EDGetTokenT<HGCalSoARecHitsExtraHostCollection> recHitsExtra;
    bool hasTiming;
  };

  edm::EDGetTokenT<reco::CaloClusterHostCollection> const getTokenSoAClusters_;
  edm::EDGetTokenT<std::vector<uint32_t>> const getTokenClusterOffsets_;
  std::vector<Input> inputs_;
  unsigned int hitsTime_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClustersFromSoAProducer);
