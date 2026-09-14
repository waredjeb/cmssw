// Counts what is actually IN the truth-branch association maps.
//
// A product being present proves only that the producer ran and put something. This
// reports rows, non-empty rows and total entries per map, which is what distinguishes
// a working associator from one that quietly wrote empty maps on every event.

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "PhysicsTools/TruthInfo/interface/SubgraphHitView.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"
#include "SimDataFormats/TruthInfo/interface/Vertex.h"

namespace {
  using SharedHitsMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  using FractionMap = ticl::TICLAssociationMap<ticl::mapWithFractionAndScore>;

  template <typename MAP>
  void report(edm::Event const& event, std::vector<std::pair<std::string, edm::EDGetTokenT<MAP>>> const& tokens) {
    for (auto const& [name, token] : tokens) {
      edm::Handle<MAP> handle;
      event.getByToken(token, handle);
      if (!handle.isValid()) {
        edm::LogPrint("TruthAssoc") << "  " << name << ": PRODUCT NOT FOUND";
        continue;
      }
      auto const& map = handle->getMap();
      std::size_t nonEmpty = 0;
      std::size_t entries = 0;
      float bestScore = std::numeric_limits<float>::infinity();
      for (auto const& row : map) {
        if (!row.empty()) {
          ++nonEmpty;
          entries += row.size();
          bestScore = std::min(bestScore, row[0].score());
        }
      }
      edm::LogPrint("TruthAssoc") << "  " << name << ": rows=" << map.size() << " nonEmpty=" << nonEmpty
                                  << " entries=" << entries
                                  << (nonEmpty > 0 ? "  bestScore=" + std::to_string(bestScore) : "");
    }
  }

  // The per-object detail behind the row counts above, for CALORIMETRIC domains whose reco
  // objects are CaloClusters (PFClusters, layer clusters): the cluster itself, then each
  // matched branch with what the matcher compared. Energies on the truth side are sim-hit
  // energies from the hit index, in the calorimeter the cluster sits in and in the whole
  // Calo channel, because that is what the shared energy and the score are made of; the
  // branch four-momentum is the generator or SimTrack momentum of the root. The hits are
  // read through SubgraphHitView, the accessor the associator's own root store agrees
  // with: the plain LogicalGraphHitIndex::subgraphHits is empty for a generator-only
  // ancestor under the shared layout, which would print zero sim energy for the very
  // roots that trivially cover every cluster.
  struct BranchSummary {
    double simEnergyInDetector = 0.;
    double simEnergyCalo = 0.;
    std::size_t cellsInDetector = 0;
  };

  BranchSummary summarizeBranch(truth::SubgraphHitView& hits, uint32_t root, DetId::Detector det) {
    BranchSummary out;
    for (auto const& hit : hits.subgraphHits(truth::HitChannel::Calo, root)) {
      out.simEnergyCalo += hit.energy;
      if (DetId(hit.detId).det() == det) {
        out.simEnergyInDetector += hit.energy;
        ++out.cellsInDetector;
      }
    }
    return out;
  }

  std::string describeRoot(truth::Graph const& graph, truth::SubgraphHitView& hits, uint32_t root, DetId::Detector det) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(2);
    if (root >= graph.nParticles()) {
      os << "root " << root << " (out of range)";
      return os.str();
    }
    auto const& data = graph.particles()[root];
    os << "root " << root << " pdg " << data.pdgId << " E=" << data.momentum.energy() << " pt=" << data.momentum.pt()
       << " eta=" << data.momentum.eta() << " phi=" << data.momentum.phi();
    const auto production = truth::Particle(&graph, root).productionVertices();
    if (!production.empty()) {
      auto const& pos = production.front().position();
      os << " vtx=(" << pos.x() << "," << pos.y() << "," << pos.z() << ")";
    }
    const auto summary = summarizeBranch(hits, root, det);
    os << " simE[det]=" << std::setprecision(4) << summary.simEnergyInDetector << " (" << summary.cellsInDetector
       << " cells) simE[calo]=" << summary.simEnergyCalo;
    return os.str();
  }

  std::string describeCluster(reco::CaloCluster const& cluster, std::size_t index) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(2);
    auto const& pos = cluster.position();
    os << "reco[" << index << "] E=" << cluster.energy() << " eta=" << pos.eta() << " phi=" << pos.phi()
       << " r=" << std::sqrt(pos.x() * pos.x() + pos.y() * pos.y()) << " z=" << pos.z()
       << " cells=" << cluster.hitsAndFractions().size() << " det=" << DetId(cluster.seed()).det();
    return os.str();
  }
}  // namespace

class TruthBranchAssociationDumper : public edm::one::EDAnalyzer<> {
public:
  explicit TruthBranchAssociationDumper(edm::ParameterSet const& cfg) {
    for (auto const& tag : cfg.getParameter<std::vector<edm::InputTag>>("sharedHitsMaps")) {
      hitsTokens_.emplace_back(tag.encode(), consumes<SharedHitsMap>(tag));
    }
    for (auto const& tag : cfg.getParameter<std::vector<edm::InputTag>>("fractionMaps")) {
      fractionTokens_.emplace_back(tag.encode(), consumes<FractionMap>(tag));
    }
    // Composite domains are only meaningful if their constituents exist: a primary
    // vertex with no tracks is correctly skipped, and that must not be mistaken for a
    // broken association.
    for (auto const& tag : cfg.getParameter<std::vector<edm::InputTag>>("vertexDiagnostics")) {
      vertexTokens_.emplace_back(tag.encode(), consumes<std::vector<reco::Vertex>>(tag));
    }
    for (auto const& pset : cfg.getParameter<std::vector<edm::ParameterSet>>("caloDetails")) {
      CaloDetail detail;
      detail.name = pset.getParameter<edm::InputTag>("clusters").encode();
      detail.clusters = consumes<edm::View<reco::CaloCluster>>(pset.getParameter<edm::InputTag>("clusters"));
      for (auto const& tag : pset.getParameter<std::vector<edm::InputTag>>("recoToTruthMaps")) {
        detail.recoToTruth.emplace_back(tag.instance(), consumes<SharedHitsMap>(tag));
      }
      detail.truthToReco = consumes<SharedHitsMap>(pset.getParameter<edm::InputTag>("truthToRecoMap"));
      caloDetails_.push_back(std::move(detail));
    }
    if (!caloDetails_.empty()) {
      graphToken_ = consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("src"));
      hitIndexToken_ = consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"));
      maxObjects_ = cfg.getParameter<unsigned int>("maxObjects");
      maxMatches_ = cfg.getParameter<unsigned int>("maxMatches");
      if (auto const tag = cfg.getParameter<edm::InputTag>("truthTargets"); !tag.label().empty()) {
        truthTargetsToken_ = consumes<std::vector<unsigned int>>(tag);
      }
    }
  }

  void analyze(edm::Event const& event, edm::EventSetup const&) override {
    edm::LogPrint("TruthAssoc") << "=== event " << event.id().event() << " ===";
    report(event, hitsTokens_);
    report(event, fractionTokens_);
    for (auto const& [name, token] : vertexTokens_) {
      edm::Handle<std::vector<reco::Vertex>> handle;
      event.getByToken(token, handle);
      if (!handle.isValid()) {
        continue;
      }
      for (std::size_t i = 0; i < handle->size(); ++i) {
        auto const& v = (*handle)[i];
        edm::LogPrint("TruthAssoc") << "  " << name << "[" << i << "]: tracks=" << v.tracksSize()
                                    << " isFake=" << v.isFake() << " ndof=" << v.ndof();
      }
    }
    if (!caloDetails_.empty()) {
      reportCaloDetails(event);
    }
  }

  void reportCaloDetails(edm::Event const& event) const {
    auto const& graph = event.get(graphToken_);
    truth::SubgraphHitView hitIndex(event.get(hitIndexToken_));
    // Optional restriction of the truth-driven print to one level's target list, so the
    // roots shown are the ones an efficiency at that level is counted over rather than
    // the generator ancestors that come first by index.
    std::vector<bool> isTarget;
    if (!truthTargetsToken_.isUninitialized()) {
      const edm::Handle<std::vector<unsigned int>> targets = event.getHandle(truthTargetsToken_);
      if (targets.isValid()) {
        isTarget.assign(graph.nParticles(), false);
        for (const unsigned int id : *targets) {
          if (id < isTarget.size()) {
            isTarget[id] = true;
          }
        }
      }
    }
    for (auto const& detail : caloDetails_) {
      const edm::Handle<edm::View<reco::CaloCluster>> clusters = event.getHandle(detail.clusters);
      if (!clusters.isValid()) {
        edm::LogPrint("TruthAssoc") << "  " << detail.name << ": CLUSTERS NOT FOUND";
        continue;
      }
      edm::LogPrint("TruthAssoc") << "  --- " << detail.name << ": reco-driven, first " << maxObjects_
                                  << " clusters, up to " << maxMatches_ << " matches each ---";
      const std::size_t nShown = std::min<std::size_t>(clusters->size(), maxObjects_);
      for (std::size_t i = 0; i < nShown; ++i) {
        auto const& cluster = (*clusters)[i];
        const auto det = DetId(cluster.seed()).det();
        edm::LogPrint("TruthAssoc") << "  " << describeCluster(cluster, i);
        for (auto const& [wp, token] : detail.recoToTruth) {
          const edm::Handle<SharedHitsMap> map = event.getHandle(token);
          if (!map.isValid() || i >= map->getMap().size()) {
            continue;
          }
          auto const& row = map->getMap()[i];
          if (row.empty()) {
            edm::LogPrint("TruthAssoc") << "      " << wp << ": no match";
            continue;
          }
          const std::size_t nMatches = std::min<std::size_t>(row.size(), maxMatches_);
          for (std::size_t m = 0; m < nMatches; ++m) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(4) << "      " << wp << " [" << m << "/" << row.size()
               << "] sharedSimE=" << row[m].sharedEnergy() << " score=" << row[m].score() << "  "
               << describeRoot(graph, hitIndex, row[m].index(), det);
            edm::LogPrint("TruthAssoc") << os.str();
          }
        }
      }
      const edm::Handle<SharedHitsMap> truthToReco = event.getHandle(detail.truthToReco);
      if (!truthToReco.isValid()) {
        continue;
      }
      edm::LogPrint("TruthAssoc") << "  --- " << detail.name << ": truth-driven, first " << maxObjects_
                                  << " matched roots" << (isTarget.empty() ? "" : " among the level targets") << " ---";
      std::size_t shown = 0;
      auto const& rows = truthToReco->getMap();
      for (std::size_t root = 0; root < rows.size() && shown < maxObjects_; ++root) {
        if (rows[root].empty() || (!isTarget.empty() && !isTarget[root])) {
          continue;
        }
        ++shown;
        // The detector for the sim-energy breakdown is the one the best-matched cluster sits in.
        const std::size_t best = rows[root][0].index();
        const auto det = best < clusters->size() ? DetId((*clusters)[best].seed()).det() : DetId::Detector(0);
        edm::LogPrint("TruthAssoc") << "  " << describeRoot(graph, hitIndex, static_cast<uint32_t>(root), det);
        const std::size_t nMatches = std::min<std::size_t>(rows[root].size(), maxMatches_);
        for (std::size_t m = 0; m < nMatches; ++m) {
          auto const& e = rows[root][m];
          std::ostringstream os;
          os << std::fixed << std::setprecision(4) << "      [" << m << "/" << rows[root].size()
             << "] sharedEnergyFraction=" << e.sharedEnergy() << " reverseScore=" << e.score() << "  ";
          if (e.index() < clusters->size()) {
            os << describeCluster((*clusters)[e.index()], e.index());
          } else {
            os << "reco[" << e.index() << "] out of range";
          }
          edm::LogPrint("TruthAssoc") << os.str();
        }
      }
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::vector<edm::InputTag>>("sharedHitsMaps", {});
    desc.add<std::vector<edm::InputTag>>("fractionMaps", {});
    desc.add<std::vector<edm::InputTag>>("vertexDiagnostics", {});
    edm::ParameterSetDescription calo;
    calo.add<edm::InputTag>("clusters")->setComment("a CaloCluster-derived collection, e.g. reco::PFCluster");
    calo.add<std::vector<edm::InputTag>>("recoToTruthMaps", {})
        ->setComment("RecoToTruth maps of that collection, one per working point to print");
    calo.add<edm::InputTag>("truthToRecoMap")->setComment("the TruthToReco map of that collection");
    desc.addVPSet("caloDetails", calo, {})
        ->setComment(
            "Per-object print of what the matcher compared, for calorimetric domains: the cluster energy, "
            "position and cell count, then each matched branch with its pdgId, four-momentum, production "
            "vertex, sim energy in the cluster's detector and in the whole Calo channel, and the map's shared "
            "energy and score. Needs src and hitIndex");
    desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"));
    desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
    desc.add<unsigned int>("maxObjects", 10)->setComment("caloDetails: clusters and matched roots shown per event");
    desc.add<unsigned int>("maxMatches", 3)->setComment("caloDetails: matches shown per object");
    desc.add<edm::InputTag>("truthTargets", edm::InputTag())
        ->setComment(
            "caloDetails: a level target list from TruthBranchTargetsProducer, e.g. "
            "truthBranchTargets:truthToRecoTargetsCaloBoundary, to restrict the truth-driven print to the roots "
            "an efficiency at that level counts. Empty label: every root with a match, ancestors first");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  std::vector<std::pair<std::string, edm::EDGetTokenT<SharedHitsMap>>> hitsTokens_;
  std::vector<std::pair<std::string, edm::EDGetTokenT<FractionMap>>> fractionTokens_;
  std::vector<std::pair<std::string, edm::EDGetTokenT<std::vector<reco::Vertex>>>> vertexTokens_;

  struct CaloDetail {
    std::string name;
    edm::EDGetTokenT<edm::View<reco::CaloCluster>> clusters;
    std::vector<std::pair<std::string, edm::EDGetTokenT<SharedHitsMap>>> recoToTruth;
    edm::EDGetTokenT<SharedHitsMap> truthToReco;
  };
  std::vector<CaloDetail> caloDetails_;
  edm::EDGetTokenT<truth::Graph> graphToken_;
  edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  edm::EDGetTokenT<std::vector<unsigned int>> truthTargetsToken_;
  unsigned int maxObjects_ = 10;
  unsigned int maxMatches_ = 3;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TruthBranchAssociationDumper);
