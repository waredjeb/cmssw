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
// It also books, through TFileService, kinematic histograms (eta, energy, pt, phi) of
// the truth particles, split by PDG category, for several populations:
//   gen            every particle with a generator record
//   genStable      the generator-stable ones (status 1)
//   <level>        the level antichain of each configured truth level, straight from the
//                  graph with no kinematic selection, e.g. reconstructableFromSignal
// Output layout: <population>/<pdgCategory>/{eta,energy,pt,phi}, plus an "all" category.
//
// Unlike a DQMGlobalEDAnalyzer this is an edm::one module: it sees one event at a time,
// so per-event scratch objects such as truth::SubgraphHitView can be plain locals.

#include <algorithm>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TH1F.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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

namespace {
  // PDG categories the kinematic histograms are split by, keyed by |pdgId|. Anything
  // else lands in "other"; every particle also fills "all".
  const std::vector<std::pair<int, const char*>> kPdgCategories = {
      {11, "electron"}, {13, "muon"},     {15, "tau"},   {22, "photon"},   {111, "pi0"},      {211, "pion"},
      {130, "K0L"},     {310, "K0S"},     {321, "kaon"}, {2212, "proton"}, {2112, "neutron"}, {12, "neutrino"},
      {14, "neutrino"}, {16, "neutrino"}, {21, "gluon"}, {1, "quark"},     {2, "quark"},      {3, "quark"},
      {4, "quark"},     {5, "quark"},     {6, "top"},    {23, "Z"},        {24, "W"},         {25, "Higgs"}};

  const char* pdgCategory(int pdgId) {
    const int absId = std::abs(pdgId);
    for (auto const& [id, name] : kPdgCategories) {
      if (id == absId) {
        return name;
      }
    }
    return "other";
  }

  struct KinematicHistograms {
    TH1F* eta = nullptr;
    TH1F* energy = nullptr;
    TH1F* pt = nullptr;
    TH1F* phi = nullptr;

    void fill(math::XYZTLorentzVectorD const& p4) const {
      eta->Fill(p4.eta());
      energy->Fill(p4.energy());
      pt->Fill(p4.pt());
      phi->Fill(p4.phi());
    }
  };
}  // namespace

class TruthGraphPFClusterExplorer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
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

  // Kinematic histograms: population -> pdg category -> the four histograms. Booked in
  // the constructor for every category, so a category that never fills stays empty
  // rather than absent and the directory layout is the same on every sample.
  std::vector<truth::Level> histogramLevels_;
  std::map<std::string, std::map<std::string, KinematicHistograms>> kinematics_;
  void bookPopulation(TFileService& fs, std::string const& population);
  void fillPopulation(std::string const& population, truth::Graph const& graph, uint32_t particleId);
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
      domain.recoToTruth.emplace_back(wp,
                                      consumes<CaloMap>(edm::InputTag(associator, domain.key + "RecoToTruth" + wp)));
    }
    domain.truthToReco = consumes<CaloMap>(edm::InputTag(associator, domain.key + "TruthToReco"));
    domains_.push_back(std::move(domain));
  }

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  bookPopulation(*fs, "gen");
  bookPopulation(*fs, "genStable");
  for (auto const& name : cfg.getParameter<std::vector<std::string>>("histogramLevels")) {
    histogramLevels_.push_back(truth::levelFromName(name));
    bookPopulation(*fs, name);
  }
}

void TruthGraphPFClusterExplorer::bookPopulation(TFileService& fs, std::string const& population) {
  std::vector<std::string> categories = {"all", "other"};
  for (auto const& [id, name] : kPdgCategories) {
    if (std::find(categories.begin(), categories.end(), name) == categories.end()) {
      categories.emplace_back(name);
    }
  }
  TFileDirectory dir = fs.mkdir(population);
  for (auto const& category : categories) {
    TFileDirectory sub = dir.mkdir(category);
    KinematicHistograms& h = kinematics_[population][category];
    const std::string title = population + " " + category;
    h.eta = sub.make<TH1F>("eta", (title + ";#eta;particles").c_str(), 100, -5., 5.);
    h.energy = sub.make<TH1F>("energy", (title + ";E [GeV];particles").c_str(), 200, 0., 500.);
    h.pt = sub.make<TH1F>("pt", (title + ";p_{T} [GeV];particles").c_str(), 200, 0., 200.);
    h.phi = sub.make<TH1F>("phi", (title + ";#phi;particles").c_str(), 64, -3.2, 3.2);
  }
}

void TruthGraphPFClusterExplorer::fillPopulation(std::string const& population,
                                                 truth::Graph const& graph,
                                                 uint32_t particleId) {
  auto const& data = graph.particles()[particleId];
  auto& byCategory = kinematics_[population];
  byCategory["all"].fill(data.momentum);
  byCategory[pdgCategory(data.pdgId)].fill(data.momentum);
}

void TruthGraphPFClusterExplorer::analyze(edm::Event const& event, edm::EventSetup const&) {
  // ---- truth side -----------------------------------------------------------------
  auto const& graph = event.get(graphToken_);

  // Kinematics of the truth populations. GEN particles are those with a generator
  // record; the synthetic ones the graph invents (interaction nodes, collapsed pileup
  // vertices) carry an accounting momentum and are skipped. A level antichain comes
  // straight from the graph, with no pt/eta selection: reconstructableFromSignal needs
  // the Signal flag stamped at DIGI and is empty on a sample produced without it, where
  // reconstructableFinalState is the event-wide equivalent.
  for (uint32_t id = 0; id < graph.nParticles(); ++id) {
    auto const& data = graph.particles()[id];
    if (!data.hasGen() || data.isSynthetic()) {
      continue;
    }
    fillPopulation("gen", graph, id);
    if (data.status == 1) {
      fillPopulation("genStable", graph, id);
    }
  }
  for (const truth::Level level : histogramLevels_) {
    for (const uint32_t id : truth::levelAntichain(graph, level)) {
      fillPopulation(truth::levelName(level), graph, id);
    }
  }
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
    auto const& truthToReco = event.get(domain.truthToReco).getMap();                           // rows = particles
    std::vector<std::pair<std::string, ticl::mapWithSharedEnergyAndScore const*>> recoToTruth;  // rows = clusters
    for (auto const& [wp, token] : domain.recoToTruth) {
      recoToTruth.emplace_back(wp, &event.get(token).getMap());
    }

    //
    // Reading the maps:
    //   for (std::size_t i = 0; i < clusters.size(); ++i){
    //     for (auto const& e : (*recoToTruth[wpIndex].second)[i]){   // best first
    //       e.index()         //-> root particle id
    //       e.sharedEnergy()  //-> shared SIM energy (sum over cells of min(reco fraction, branch) x cell energy)
    //       e.score()         //-> reco-normalised score, 0 = branch covers the whole cluster
    //      }
    //    }
    //   for (unsigned int root : *levelTargets[k].second)
    //     for (auto const& e : truthToReco[root])
    //       e.index()         //-> cluster index
    //       e.sharedEnergy()  //-> shared energy FRACTION of the branch energy in the denominator detectors
    //       e.score()         //-> branch-normalised (reverse) score over the whole Calo channel

    // Fixed lists every sharing root; the adaptive working points keep one per cluster.

    // ---- your study goes here ---------------------------------------------------------
    edm::LogInfo("TruthGraphPFClusterExplorer")
        << domain.key << ": " << clusters.size() << " clusters, " << recoToTruth.size() << " reco-driven maps, "
        << truthToReco.size() << " truth rows, " << selectedRoots.size() << " candidate roots, " << levelTargets.size()
        << " levels, " << hits.nParticles() << " particles";
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
  desc.add<std::vector<std::string>>("histogramLevels",
                                     {"reconstructableFromSignal", "reconstructableFinalState", "caloBoundary"})
      ->setComment("truth levels whose antichain gets its own kinematic histogram population");
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
