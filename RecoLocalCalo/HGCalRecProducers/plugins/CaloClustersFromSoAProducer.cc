// Rebuild the legacy AoS layer clusters from the portable SoA representation.
//
// This is the single compatibility bridge between the SoA world (layer-cluster
// production, merging and TICL) and the consumers that need `edm::Ref`/`edm::Ptr`
// into an addressable collection and therefore cannot read an SoA at all:
// the sim-truth associators, HGCal validation, PF and the EGamma HLT producers.
//
// It sits once, after the merge, and reconstructs the hits and fractions from
// the association map, so it works for every subdetector that feeds the merge
// (EE, FH, BH and, under `ticl_barrel`, EB/HB) rather than only for the ones
// whose rechit SoA happens to be available.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"

#include <memory>
#include <utility>
#include <vector>

class CaloClustersFromSoAProducer : public edm::stream::EDProducer<> {
public:
  CaloClustersFromSoAProducer(edm::ParameterSet const& config)
      // The cluster SoA and its association map are emitted by the same module
      // under the same (unnamed) instance: EDM resolves them by type.
      : clustersToken_(consumes(config.getParameter<edm::InputTag>("src"))),
        hitsAndFractionsToken_(consumes(config.getParameter<edm::InputTag>("src"))),
        timeClname_(config.getParameter<std::string>("timeClname")) {
    produces<std::vector<reco::BasicCluster>>();
    produces<edm::ValueMap<std::pair<float, float>>>(timeClname_);
  }

  ~CaloClustersFromSoAProducer() override = default;

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
    auto const clusters_v = iEvent.get(clustersToken_).const_view();
    auto const hits_v = iEvent.get(hitsAndFractionsToken_).const_view();
    const int numberOfClusters = clusters_v.position().metadata().size();

    auto clusters = std::make_unique<std::vector<reco::BasicCluster>>();
    clusters->reserve(numberOfClusters);
    std::vector<std::pair<float, float>> times;
    times.reserve(numberOfClusters);

    for (int i = 0; i < numberOfClusters; ++i) {
      auto const hitsOfCluster = hits_v[i];
      std::vector<std::pair<DetId, float>> hitsAndFractions;
      hitsAndFractions.reserve(hitsOfCluster.size());
      for (auto const& hitAndFraction : hitsOfCluster) {
        hitsAndFractions.emplace_back(hitAndFraction.hit, hitAndFraction.fraction);
      }

      clusters->emplace_back(
          clusters_v.energy()[i].energy(),
          math::XYZPoint(clusters_v.position()[i].x(), clusters_v.position()[i].y(), clusters_v.position()[i].z()),
          clusters_v.indexes()[i].caloID(),
          hitsAndFractions,
          clusters_v.indexes()[i].algoID(),
          clusters_v.indexes()[i].seedID(),
          clusters_v.indexes()[i].flags());
      // The reco::CaloCluster constructor defaults both to -1; copy them across
      // anyway so the converter stays faithful if a producer ever sets them.
      clusters->back().setCorrectedEnergy(clusters_v.energy()[i].correctedEnergy());
      clusters->back().setCorrectedEnergyUncertainty(clusters_v.energy()[i].correctedEnergyUncertainty());

      times.emplace_back(clusters_v.timing()[i].time(), clusters_v.timing()[i].timeError());
    }

    auto clusterHandle = iEvent.put(std::move(clusters));

    auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
    edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
    filler.insert(clusterHandle, times.begin(), times.end());
    filler.fill();
    iEvent.put(std::move(timeCl), timeClname_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hgcalMergeLayerClusters"))
        ->setComment("module emitting both the cluster SoA and the hits-and-fractions association map");
    desc.add<std::string>("timeClname", "timeLayerCluster");
    descriptions.add("caloClustersFromSoA", desc);
  }

private:
  edm::EDGetTokenT<reco::CaloClusterHostCollection> const clustersToken_;
  edm::EDGetTokenT<ticl::HitsAndFractionsHost> const hitsAndFractionsToken_;
  std::string const timeClname_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloClustersFromSoAProducer);
