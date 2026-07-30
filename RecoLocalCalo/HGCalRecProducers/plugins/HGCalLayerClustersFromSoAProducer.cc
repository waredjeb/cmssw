#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsExtraHostCollection.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"

#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include <vector>

class HGCalLayerClustersFromSoAProducer : public edm::stream::EDProducer<> {
public:
  HGCalLayerClustersFromSoAProducer(edm::ParameterSet const& config)
      : getTokenSoAClusters_(consumes(config.getParameter<edm::InputTag>("src"))),
        getTokenSoACells_(consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))),
        getTokenSoARecHitsExtra_(consumes(config.getParameter<edm::InputTag>("hgcalRecHitsLayerClustersSoA"))),
        detector_(config.getParameter<std::string>("detector")),
        hitsTime_(config.getParameter<unsigned int>("nHitsTime")) {
    if (detector_ == "HFNose") {
      algoId_ = reco::CaloCluster::hfnose;
    } else if (detector_ == "EE") {
      algoId_ = reco::CaloCluster::hgcal_em;
    } else if (detector_ == "BH") {
      algoId_ = reco::CaloCluster::hgcal_scintillator;
    } else {  //for FH
      algoId_ = reco::CaloCluster::hgcal_had;
    }

    produces<std::vector<float>>("InitialLayerClustersMask");

    produces<reco::CaloClusterHostCollection>();
    produces<ticl::HitsAndFractionsHost>();
  }

  ~HGCalLayerClustersFromSoAProducer() override = default;

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
    auto const& deviceData = iEvent.get(getTokenSoAClusters_);

    auto const& deviceSoARecHitsExtra = iEvent.get(getTokenSoARecHitsExtra_);
    auto const soaRecHitsExtra_v = deviceSoARecHitsExtra.view();

    auto const& deviceSoACells = iEvent.get(getTokenSoACells_);
    auto const soaCells_v = deviceSoACells.view();

    auto const deviceView = deviceData.view();
    auto const position_v = deviceView.position();
    auto const energy_v = deviceView.energy();
    auto const indexes_v = deviceView.indexes();
    const int numberOfClusters = position_v.metadata().size();
    const int numberOfRecHits = soaRecHitsExtra_v.metadata().size();

    // Number of rechits in each cluster, excluding outliers
    std::vector<int> hitsPerCluster(numberOfClusters, 0);
    int numberOfClusteredHits = 0;
    for (int i = 0; i < numberOfRecHits; ++i) {
      const auto clusterIndex = soaRecHitsExtra_v[i].clusterIndex();
      if (clusterIndex == -1) {
        continue;
      }
      ++hitsPerCluster[clusterIndex];
      ++numberOfClusteredHits;
    }

    auto clusters = std::make_unique<reco::CaloClusterHostCollection>(
        cms::alpakatools::host(), numberOfClusters, numberOfClusters, numberOfClusters, numberOfClusters);
    auto hitsAndFractions = std::make_unique<ticl::HitsAndFractionsHost>(
        cms::alpakatools::host(), numberOfClusteredHits, numberOfClusters);
    auto clusters_v = clusters->view();
    auto hits_v = hitsAndFractions->view();

    // Prefix-sum the per-cluster counts into the CSR offsets, and copy the
    // per-cluster features across. algoID/caloID are reqritten from this
    // module's own detector setting rather than taken from the kernel, which
    // hardcodes the EE values.
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
      clusters_v.indexes()[i].algoID() = algoId_;
      clusters_v.indexes()[i].seedID() = indexes_v[i].seedID();
      clusters_v.indexes()[i].flags() = indexes_v[i].flags();
    }
    hits_v.offsets()[numberOfClusters].keys_offsets() = hitsOffset;
    assert(hitsOffset == numberOfClusteredHits);

    // scatter the hits, and collect the per-hit times needed below.
    std::vector<int> cursor(numberOfClusters, 0);
    std::vector<std::vector<float>> times(numberOfClusters);
    std::vector<std::vector<float>> timeErrors(numberOfClusters);
    for (int i = 0; i < numberOfRecHits; ++i) {
      const auto clusterIndex = soaRecHitsExtra_v[i].clusterIndex();
      if (clusterIndex == -1) {
        continue;
      }
      const auto slot = hits_v.offsets()[clusterIndex].keys_offsets() + cursor[clusterIndex]++;
      hits_v.content().values()[slot] = ticl::HitAndFraction{soaCells_v[i].detid(), 1.f};
      if (soaCells_v[i].timeError() < 0.f) {
        continue;
      }
      times[clusterIndex].push_back(soaCells_v[i].time());
      timeErrors[clusterIndex].push_back(1.f / (soaCells_v[i].timeError() * soaCells_v[i].timeError()));
    }

    // Assign time to each cluster
    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    for (int i = 0; i < numberOfClusters; ++i) {
      const auto timeCl = (detector_ != "BH")
                              ? timeEstimator.fixSizeHighestDensity(times[i], timeErrors[i], hitsTime_)
                              : std::pair<float, float>(-99.f, -1.f);
      clusters_v.timing()[i].time() = timeCl.first;
      clusters_v.timing()[i].timeError() = timeCl.second;
    }

    // HFNoise has no mergint step, so create the layerClustersMask here
    if (detector_ == "HFNose") {
      std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
      layerClustersMask->resize(numberOfClusters, 1.0);
      iEvent.put(std::move(layerClustersMask), "InitialLayerClustersMask");
    }

    iEvent.put(std::move(clusters));
    iEvent.put(std::move(hitsAndFractions));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hltHgcalSoALayerClustersProducer"));
    desc.add<edm::InputTag>("hgcalRecHitsLayerClustersSoA", edm::InputTag("hltHgcalSoARecHitsLayerClustersProducer"));
    desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("hltHgcalSoARecHitsProducer"));
    desc.add<unsigned int>("nHitsTime", 3);
    desc.add<std::string>("timeClname", "timeLayerCluster");
    desc.ifValue(edm::ParameterDescription<std::string>(
                     "detector", "EE", true, edm::Comment("the HGCAL component used to create clusters.")),
                 edm::allowedValues<std::string>("EE", "FH", "BH"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  edm::EDGetTokenT<reco::CaloClusterHostCollection> const getTokenSoAClusters_;
  edm::EDGetTokenT<HGCalSoARecHitsHostCollection> const getTokenSoACells_;
  edm::EDGetTokenT<HGCalSoARecHitsExtraHostCollection> const getTokenSoARecHitsExtra_;
  std::string detector_;
  unsigned int hitsTime_;
  reco::CaloCluster::AlgoId algoId_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClustersFromSoAProducer);
