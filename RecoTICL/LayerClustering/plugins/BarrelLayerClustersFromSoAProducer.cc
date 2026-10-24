#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

// PFRecHitHostCollection pulls in Eigen, which must precede any other SoA header
#include "DataFormats/ParticleFlowReco/interface/PFRecHitHostCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsExtraHostCollection.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"
#include "RecoParticleFlow/PFClusterProducer/interface/CaloRecHitResolutionProvider.h"

class BarrelLayerClustersFromSoAProducer : public edm::stream::EDProducer<> {
public:
  BarrelLayerClustersFromSoAProducer(edm::ParameterSet const& config)
      : getTokenSoAClusters_(consumes(config.getParameter<edm::InputTag>("src"))),
        getTokenSoARecHits_(consumes(config.getParameter<edm::InputTag>("recHits"))),
        getTokenSoARecHitsExtra_(consumes(config.getParameter<edm::InputTag>("barrelRecHitsLayerClustersSoA"))),
        nHitsTime_(config.getParameter<unsigned int>("nHitsTime")),
        timeClname_(config.getParameter<std::string>("timeClname")),
        timeResolutionCalc_(config.getParameterSet("timeResolutionCalc")) {
    if (config.getParameter<std::string>("detector") == "EB") {
      caloId_ = reco::CaloID::DET_ECAL_BARREL;
      algoId_ = reco::CaloCluster::barrel_em;
    } else {
      caloId_ = reco::CaloID::DET_HCAL_BARREL;
      algoId_ = reco::CaloCluster::barrel_had;
    }

    produces<std::vector<float>>("InitialLayerClustersMask");
    produces<std::vector<reco::BasicCluster>>();
    produces<edm::ValueMap<std::pair<float, float>>>(timeClname_);
  }

  ~BarrelLayerClustersFromSoAProducer() override = default;

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
    auto const& soaClusters = iEvent.get(getTokenSoAClusters_);
    auto const clusters_v = soaClusters.const_view();
    auto const& soaRecHits = iEvent.get(getTokenSoARecHits_);
    auto const recHits_v = soaRecHits.const_view();
    auto const& soaRecHitsExtra = iEvent.get(getTokenSoARecHitsExtra_);
    auto const recHitsExtra_v = soaRecHitsExtra.const_view();

    const auto nClusters = soaClusters->metadata().size()[0u];
    auto clusters = std::make_unique<std::vector<reco::BasicCluster>>();
    clusters->reserve(nClusters);

    // Times and inverse squared time resolutions of the rechits of each cluster
    std::vector<std::vector<float>> times(nClusters);
    std::vector<std::vector<float>> timeErrors(nClusters);

    for (auto i = 0; i < nClusters; ++i) {
      std::vector<std::pair<DetId, float>> thisCluster;
      thisCluster.reserve(clusters_v.position()[i].cells());
      clusters->emplace_back(
          clusters_v.energy()[i].energy(),
          math::XYZPoint(clusters_v.position()[i].x(), clusters_v.position()[i].y(), clusters_v.position()[i].z()),
          caloId_,
          std::move(thisCluster),
          algoId_);
      clusters->back().setSeed(clusters_v.indexes()[i].seedID());
      times[i].reserve(clusters_v.position()[i].cells());
      timeErrors[i].reserve(clusters_v.position()[i].cells());
    }

    // The per-rechit clustering output has the same cardinality as the PFRecHits
    for (auto i = 0; i < recHitsExtra_v.metadata().size(); ++i) {
      const int clusterIndex = recHitsExtra_v[i].clusterIndex();
      if (clusterIndex == -1) {
        continue;
      }
      assert(clusterIndex < static_cast<int>(clusters->size()));
      // Without sharing every rechit belongs to a single cluster
      (*clusters)[clusterIndex].addHitAndFraction(recHits_v[i].detId(), 1.f);
      const float rhTimeE = timeResolutionCalc_.timeResolution2(recHits_v[i].energy());
      if (rhTimeE < 0.f) {
        continue;
      }
      times[clusterIndex].push_back(recHits_v[i].time());
      timeErrors[clusterIndex].push_back(1.f / rhTimeE);
    }

    // As in the legacy BarrelLayerClusterProducer, the time is only computed for clusters with enough rechits
    std::vector<std::pair<float, float>> clusterTimes;
    clusterTimes.reserve(clusters->size());
    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    for (unsigned i = 0; i < clusters->size(); ++i) {
      if ((*clusters)[i].size() >= nHitsTime_) {
        clusterTimes.push_back(timeEstimator.fixSizeHighestDensity(times[i], timeErrors[i], nHitsTime_));
      } else {
        clusterTimes.emplace_back(-99.f, -1.f);
      }
    }

    auto clusterHandle = iEvent.put(std::move(clusters));

    auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
    edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
    filler.insert(clusterHandle, clusterTimes.begin(), clusterTimes.end());
    filler.fill();
    iEvent.put(std::move(timeCl), timeClname_);

    auto layerClustersMask = std::make_unique<std::vector<float>>(clusterHandle->size(), 1.f);
    iEvent.put(std::move(layerClustersMask), "InitialLayerClustersMask");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("TO BE DEFINED"));
    desc.add<edm::InputTag>("recHits", edm::InputTag("TO BE DEFINED"));
    desc.add<edm::InputTag>("barrelRecHitsLayerClustersSoA", edm::InputTag("TO BE DEFINED"));
    desc.ifValue(edm::ParameterDescription<std::string>("detector", "EB", true), edm::allowedValues<std::string>("EB", "HB"));

    edm::ParameterSetDescription timeResolutionCalcDesc;
    timeResolutionCalcDesc.addNode(edm::ParameterDescription<double>("noiseTerm", 1.10889, true) and
                                   edm::ParameterDescription<double>("constantTerm", 0.428192, true) and
                                   edm::ParameterDescription<double>("corrTermLowE", 0.0510871, true) and
                                   edm::ParameterDescription<double>("threshLowE", 0.5, true) and
                                   edm::ParameterDescription<double>("constantTermLowE", 0.0, true) and
                                   edm::ParameterDescription<double>("noiseTermLowE", 1.31883, true) and
                                   edm::ParameterDescription<double>("threshHighE", 5.0, true));
    desc.add<edm::ParameterSetDescription>("timeResolutionCalc", timeResolutionCalcDesc);

    desc.add<unsigned int>("nHitsTime", 3);
    desc.add<std::string>("timeClname", "timeLayerCluster");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  edm::EDGetTokenT<reco::CaloClusterHostCollection> const getTokenSoAClusters_;
  edm::EDGetTokenT<reco::PFRecHitHostCollection> const getTokenSoARecHits_;
  edm::EDGetTokenT<HGCalSoARecHitsExtraHostCollection> const getTokenSoARecHitsExtra_;
  reco::CaloID caloId_;
  reco::CaloCluster::AlgoId algoId_;
  const unsigned int nHitsTime_;
  const std::string timeClname_;
  CaloRecHitResolutionProvider timeResolutionCalc_;
};

DEFINE_FWK_MODULE(BarrelLayerClustersFromSoAProducer);
