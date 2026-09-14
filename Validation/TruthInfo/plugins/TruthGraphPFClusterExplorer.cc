// Skeleton analyzer for exploring the truth graph against particle-flow clusters.
//
// It consumes everything the PFCluster association chain produces, retrieves it in
// analyze() and stops there: the study itself is yours to write. The comments in
// analyze() show how each product is read, so the pieces can be combined without
// hunting through the producers.
//
// Products consumed:
//   - truth::Graph                    the logical truth graph (particles, vertices)
//   - truth::LogicalGraphHitIndex     per-particle direct/subgraph detector hits
//   - selectedRoots                   the branch candidates every associator matched against
//   - one target list per truth level (the efficiency denominators), by level name
//   - one or more PFCluster collections
//   - per collection: one RecoToTruth map per working point and the TruthToReco map
//
// Unlike a DQMGlobalEDAnalyzer this is an edm::one module: it sees one event at a time,
// so per-event scratch objects such as truth::SubgraphHitView can be plain locals.

#include <string>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"
#include "SimDataFormats/TruthInfo/interface/Vertex.h"

#include "PhysicsTools/TruthInfo/interface/Branch.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "PhysicsTools/TruthInfo/interface/SubgraphHitView.h"
#include "PhysicsTools/TruthInfo/interface/TruthLevels.h"

class TruthGraphPFClusterExplorer : public edm::one::EDAnalyzer<> {
public:
  explicit TruthGraphPFClusterExplorer(edm::ParameterSet const& cfg);
  void analyze(edm::Event const& event, edm::EventSetup const& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // The calorimetric association maps: value = shared energy (RecoToTruth) or shared
  // energy fraction (TruthToReco), score = reco-normalised (RecoToTruth) or
  // branch-normalised (TruthToReco). Lower score is better; row [0] is the best match.
  using CaloMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;

  struct ClusterDomain {
    std::string key;  // the collection label, also the association product key
    edm::EDGetTokenT<std::vector<reco::PFCluster>> clusters;
    // (working point name, token) in the configured order
    std::vector<std::pair<std::string, edm::EDGetTokenT<CaloMap>>> recoToTruth;
    edm::EDGetTokenT<CaloMap> truthToReco;
  };

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> selectedRootsToken_;
  // (level name, token) for the per-level efficiency denominators
  std::vector<std::pair<std::string, edm::EDGetTokenT<std::vector<unsigned int>>>> levelTargetTokens_;
  std::vector<ClusterDomain> domains_;
};

TruthGraphPFClusterExplorer::TruthGraphPFClusterExplorer(edm::ParameterSet const& cfg)
    : graphToken_(consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("src"))),
      hitIndexToken_(consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"))),
      selectedRootsToken_(consumes<std::vector<unsigned int>>(cfg.getParameter<edm::InputTag>("selectedRoots"))) {
  // One denominator per level, produced by TruthBranchTargetsProducer under the instance
  // "truthToRecoTargets" + Level (first letter capitalised), e.g. truthToRecoTargetsCaloBoundary.
  const auto targetsLabel = cfg.getParameter<std::string>("targetsProducer");
  for (auto const& level : cfg.getParameter<std::vector<std::string>>("truthLevels")) {
    std::string instance = level;
    instance[0] = std::toupper(static_cast<unsigned char>(instance[0]));
    levelTargetTokens_.emplace_back(
        level, consumes<std::vector<unsigned int>>(edm::InputTag(targetsLabel, "truthToRecoTargets" + instance)));
  }

  // Association products are keyed <collection label>[_<instance>] + RecoToTruth<WP> and
  // + TruthToReco, the same rule AllRecoToTruthBranchAssociatorsProducer uses to name them.
  const auto workingPoints = cfg.getParameter<std::vector<std::string>>("workingPoints");
  for (auto const& pset : cfg.getParameter<std::vector<edm::ParameterSet>>("domains")) {
    ClusterDomain domain;
    const auto clustersTag = pset.getParameter<edm::InputTag>("clusters");
    const auto associator = pset.getParameter<std::string>("associator");
    domain.key = clustersTag.label() + (clustersTag.instance().empty() ? "" : "_" + clustersTag.instance());
    domain.clusters = consumes<std::vector<reco::PFCluster>>(clustersTag);
    for (auto const& wp : workingPoints) {
      domain.recoToTruth.emplace_back(wp, consumes<CaloMap>(edm::InputTag(associator, domain.key + "RecoToTruth" + wp)));
    }
    domain.truthToReco = consumes<CaloMap>(edm::InputTag(associator, domain.key + "TruthToReco"));
    domains_.push_back(std::move(domain));
  }
}

void TruthGraphPFClusterExplorer::analyze(edm::Event const& event, edm::EventSetup const&) {
  // ---- truth side -----------------------------------------------------------------
  [[maybe_unused]] auto const& graph = event.get(graphToken_);
  auto const& hitIndex = event.get(hitIndexToken_);
  // Per-event coalescing view over the hit index: the accessor that is correct for EVERY
  // particle, including generator-only ancestors whose subgraph spans several ranges.
  //   hits.subgraphHits(truth::HitChannel::Calo, particleId) -> span of {detId, recHitIndex, energy, ...}
  truth::SubgraphHitView hits(hitIndex);
  // The branch candidates the associators matched against (selector-passing particles
  // at every depth, ancestors included) and the per-level efficiency denominators.
  auto const& selectedRoots = event.get(selectedRootsToken_);
  std::vector<std::pair<std::string, std::vector<unsigned int> const*>> levelTargets;
  for (auto const& [level, token] : levelTargetTokens_) {
    levelTargets.emplace_back(level, &event.get(token));
  }
  //
  // Navigating a particle:
  //   auto const& data = graph.particles()[id];          // pdgId, momentum, eventId (0 = signal)
  //   truth::Particle particle(&graph, id);               // parents(), descendants(), productionVertices()
  //   truth::Branch branch(&graph, id);                    // p4(), visibleEnergy(), stableLeaves(), isSignal()
  //   auto const& vertex = production.front().position(); // XYZT
  //
  // Recomputing an association yourself, exactly as the producer does:
  //   truth::BranchHitAssociator associator(hitIndex, selectedRoots,
  //                                         truth::BranchHitAssociator::Metric::SharedEnergy,
  //                                         truth::HitChannel::Calo, /*emptyRootsMeansAll=*/false,
  //                                         1u << DetId::Ecal);            // denominator detectors
  //   const auto recoHits = truth::recoHits(cluster);                       // (detId, fraction) cells
  //   const auto matches = associator.bestBranches(recoHits);               // all sharing roots, best first
  //   const auto adaptive = truth::BranchHitAssociator::bestAdaptiveBranch(matches, 1.f, 1.f);
  //   // BranchMatch: rootParticleId, sharedEnergy, score, reverseScore, sharedEnergyFraction

  // ---- reco side and the association maps -------------------------------------------
  for (auto const& domain : domains_) {
    auto const& clusters = event.get(domain.clusters);
    auto const& truthToReco = event.get(domain.truthToReco).getMap();  // rows = particles
    std::vector<std::pair<std::string, ticl::mapWithSharedEnergyAndScore const*>> recoToTruth;  // rows = clusters
    for (auto const& [wp, token] : domain.recoToTruth) {
      recoToTruth.emplace_back(wp, &event.get(token).getMap());
    }
    //
    // Reading the maps:
    //   for (std::size_t i = 0; i < clusters.size(); ++i)
    //     for (auto const& e : (*recoToTruth[wpIndex].second)[i])   // best first
    //       e.index()         -> root particle id
    //       e.sharedEnergy()  -> shared SIM energy (sum over cells of min(reco fraction, branch) x cell energy)
    //       e.score()         -> reco-normalised score, 0 = branch covers the whole cluster
    //   for (unsigned int root : *levelTargets[k].second)
    //     for (auto const& e : truthToReco[root])
    //       e.index()         -> cluster index
    //       e.sharedEnergy()  -> shared energy FRACTION of the branch energy in the denominator detectors
    //       e.score()         -> branch-normalised (reverse) score over the whole Calo channel
    //
    // Fixed lists every sharing root; the adaptive working points keep one per cluster.

    // ---- your study goes here ---------------------------------------------------------
    edm::LogInfo("TruthGraphPFClusterExplorer")
        << domain.key << ": " << clusters.size() << " clusters, " << recoToTruth.size() << " reco-driven maps, "
        << truthToReco.size() << " truth rows, " << selectedRoots.size() << " candidate roots, "
        << levelTargets.size() << " levels, " << hits.nParticles() << " particles";
  }
}

void TruthGraphPFClusterExplorer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<edm::InputTag>("selectedRoots", edm::InputTag("truthBranchTargets", "selectedRoots"));
  desc.add<std::string>("targetsProducer", "truthBranchTargets");
  desc.add<std::vector<std::string>>("truthLevels", {"caloBoundary", "stableDecayProducts"})
      ->setComment("levels whose truthToRecoTargets<Level> lists are read");
  desc.add<std::vector<std::string>>("workingPoints", {"Fixed", "AdaptiveTight", "AdaptiveNominal"});
  edm::ParameterSetDescription domain;
  domain.add<edm::InputTag>("clusters")->setComment("a reco::PFCluster collection");
  domain.add<std::string>("associator")->setComment("the module that produced its association maps");
  desc.addVPSet("domains",
                domain,
                {[] {
                   edm::ParameterSet p;
                   p.addParameter<edm::InputTag>("clusters", edm::InputTag("particleFlowClusterECAL"));
                   p.addParameter<std::string>("associator", "truthBranchPFClusterEcalAssociators");
                   return p;
                 }(),
                 [] {
                   edm::ParameterSet p;
                   p.addParameter<edm::InputTag>("clusters", edm::InputTag("particleFlowClusterHCAL"));
                   p.addParameter<std::string>("associator", "truthBranchPFClusterHcalAssociators");
                   return p;
                 }()});
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TruthGraphPFClusterExplorer);
