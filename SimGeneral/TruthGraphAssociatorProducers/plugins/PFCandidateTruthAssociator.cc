// Truth association and reconstruction ladder for particle-flow candidates.
//
// A PFCandidate is a composite object: a track plus calorimeter clusters in the barrel
// (through the PF blocks), a track plus tracksters in the endcap (through the pfTICL
// TICLCandidate, which carries no blocks). It is never matched on hits. Its truth is
// aggregated from the association maps its constituents already have, so a failure can
// be attributed to the step that caused it: clustering, tracking, or linking.
//
// Reco-driven (per candidate), the TRACK-FIRST rule:
// (the constituent maps only know the branch selector's candidates, so a constituent of a
// sub-threshold particle has no branch here; see PFCandidateClass::Unmatched)
//   1. a track resolved to a branch above the purity floor anchors the candidate to that
//      particle P. The calorimeter constituents decide the LEVEL: energy of P and its
//      descendants keeps the level at P; a share above maxForeignFraction belonging to
//      branches with a common ancestor A that the candidate's constituents cover well
//      lifts the level to A (merged); any other share is contamination and the level
//      stays at P.
//   2. an unmatched track falls through to 3, flagged.
//   3. without a track the level is the energy-weighted best branch of the calorimeter
//      constituents, with the same merge test between the two leading branches.
//   4. no constituent with a branch: unmatched (no constituent in any configured domain:
//      unevaluated).
//
// Truth-driven (per particle of the validation level), the LADDER: each rung stored as a
// flag with its underlying fraction, unconditionally, so a consumer draws the cut-flow
// or the per-rung efficiency. Rungs whose detector is not expected are not evaluated.
//
// Products, for the candidate collection with key K (label[_instance]):
//   K + "RecoToTruth" + wp   TICLAssociationMap, one per working point: "Fixed" holds the
//                            whole composition, any other name the resolved level
//   K + "TruthToReco"        TICLAssociationMap over the graph's particles
//   K + "TruthRecords"       vector<truth::PFCandidateTruthRecord>
//   K + "CandidateRecords"   vector<truth::PFCandidateRecord>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "HepPDT/ParticleID.hh"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/PFCandidateTruthRecords.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"

#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/SubgraphHitView.h"

namespace {
  using truth::byAscendingScore;
  using CaloMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  using Rows = ticl::mapWithSharedEnergyAndScore;
  using truth::PFCandidateClass;
  using truth::PFRung;

  constexpr uint32_t kNoBranch = truth::PFCandidateRecord::kNoBranch;

  enum class Detector : uint8_t { Ecal = 0, Hcal = 1, Hgcal = 2, kCount = 3 };
  constexpr std::size_t kDetectors = static_cast<std::size_t>(Detector::kCount);

  // The calorimeter constituent domains: one PFCluster or trackster collection each,
  // with its RecoToTruth rows at the two working points this producer reads and its
  // TruthToReco rows. Filled per event from the handles.
  struct CaloDomain {
    Detector detector;
    std::vector<float> energies;  // per constituent index
    Rows const* fixed = nullptr;  // every sharing branch, best first
    Rows const* adaptive = nullptr;
    Rows const* truthToReco = nullptr;
    edm::ProductID productId;
    bool valid() const { return fixed != nullptr && adaptive != nullptr && truthToReco != nullptr; }
  };

  // A candidate's constituents as (domain, index, energy) plus its track.
  struct Constituents {
    int32_t track = -1;
    std::vector<std::pair<std::size_t, uint32_t>> calo;  // (domain slot, constituent index)
    bool endcapMapped = false;
  };

  [[nodiscard]] Detector detectorOf(uint32_t detId) {
    switch (DetId(detId).det()) {
      case DetId::Ecal:
        return Detector::Ecal;
      case DetId::Hcal:
        return Detector::Hcal;
      default:
        return Detector::Hgcal;
    }
  }

  [[nodiscard]] bool isCharged(int pdgId) { return HepPDT::ParticleID(pdgId).threeCharge() != 0; }

  // The species table: which PF type a species should come out as. -1: not evaluated.
  [[nodiscard]] int expectedPFType(int pdgId) {
    switch (std::abs(pdgId)) {
      case 211:
      case 321:
      case 2212:
        return reco::PFCandidate::h;
      case 11:
        return reco::PFCandidate::e;
      case 13:
        return reco::PFCandidate::mu;
      case 22:
      case 111:
        return reco::PFCandidate::gamma;
      case 2112:
      case 130:
        return reco::PFCandidate::h0;
      default:
        return -1;
    }
  }
}  // namespace

class PFCandidateTruthAssociator : public edm::global::EDProducer<> {
public:
  explicit PFCandidateTruthAssociator(edm::ParameterSet const&);
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  struct CaloDomainConfig {
    Detector detector;
    edm::EDGetTokenT<edm::View<reco::CaloCluster>> collection;  // PFClusters or tracksters' proxy
    edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters;  // endcap only
    edm::EDGetTokenT<CaloMap> fixed;
    edm::EDGetTokenT<CaloMap> adaptive;
    edm::EDGetTokenT<CaloMap> truthToReco;
    bool isTrackster = false;
    double simEnergyScale = 1.;
  };

  [[nodiscard]] truth::CaloRegion regionOf(double absEta) const;
  [[nodiscard]] Constituents constituentsOf(reco::PFCandidate const& candidate,
                                            edm::ProductID const& trackProductId,
                                            std::vector<CaloDomain> const& domains,
                                            std::vector<TICLCandidate> const* ticlCandidates,
                                            std::vector<int32_t> const& ticlIndexOf,
                                            std::size_t candidateIndex) const;

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> targetsToken_;
  const edm::EDGetTokenT<std::vector<reco::PFCandidate>> candidatesToken_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<CaloMap> trackFixedToken_;
  edm::EDGetTokenT<CaloMap> trackAdaptiveToken_;
  edm::EDGetTokenT<std::vector<TICLCandidate>> ticlCandidatesToken_;
  std::vector<CaloDomainConfig> caloDomains_;
  std::string key_;
  std::vector<std::string> workingPoints_;
  std::string adaptiveName_;

  const double minTrackPurity_;
  const double minCollectedFraction_;
  const double minLinkedFraction_;
  const double maxForeignFraction_;
  const double minExpectedDetectorFraction_;
  const double minMergedCoverage_;
  // Sim-to-reco energy scale per detector, so the expected-detector rule compares
  // deposits on one scale. The hit index stores PCaloHit sampling energies: about the
  // reco energy in ECAL and HGCAL, about 1% of it in HCAL.
  std::array<double, kDetectors> simEnergyScale_{{1., 1., 1.}};
  const double barrelEtaMax_;
  const double endcapEtaMin_;
  const double endcapEtaMax_;

  mutable std::once_flag endcapWarned_;
};

PFCandidateTruthAssociator::PFCandidateTruthAssociator(edm::ParameterSet const& cfg)
    : graphToken_(consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("src"))),
      hitIndexToken_(consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"))),
      targetsToken_(consumes<std::vector<unsigned int>>(cfg.getParameter<edm::InputTag>("targets"))),
      candidatesToken_(consumes<std::vector<reco::PFCandidate>>(cfg.getParameter<edm::InputTag>("candidates"))),
      tracksToken_(consumes<std::vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracks"))),
      workingPoints_(cfg.getParameter<std::vector<std::string>>("workingPointNames")),
      adaptiveName_(cfg.getParameter<std::string>("constituentWorkingPoint")),
      minTrackPurity_(cfg.getParameter<double>("minTrackPurity")),
      minCollectedFraction_(cfg.getParameter<double>("minCollectedFraction")),
      minLinkedFraction_(cfg.getParameter<double>("minLinkedFraction")),
      maxForeignFraction_(cfg.getParameter<double>("maxForeignFraction")),
      minExpectedDetectorFraction_(cfg.getParameter<double>("minExpectedDetectorFraction")),
      minMergedCoverage_(cfg.getParameter<double>("minMergedCoverage")),
      barrelEtaMax_(cfg.getParameter<double>("barrelEtaMax")),
      endcapEtaMin_(cfg.getParameter<double>("endcapEtaMin")),
      endcapEtaMax_(cfg.getParameter<double>("endcapEtaMax")) {
  const auto candidatesTag = cfg.getParameter<edm::InputTag>("candidates");
  key_ = candidatesTag.label() + (candidatesTag.instance().empty() ? "" : "_" + candidatesTag.instance());

  // Constituent maps follow the key rule of AllRecoToTruthBranchAssociatorsProducer:
  // <collection key> + "RecoToTruth" + wp, and + "TruthToReco".
  auto mapTag = [](std::string const& module, edm::InputTag const& collection, std::string const& suffix) {
    const std::string key = collection.label() + (collection.instance().empty() ? "" : "_" + collection.instance());
    return edm::InputTag(module, key + suffix);
  };
  const auto tracksTag = cfg.getParameter<edm::InputTag>("tracks");
  const auto trackAssociator = cfg.getParameter<std::string>("trackAssociator");
  trackFixedToken_ = consumes<CaloMap>(mapTag(trackAssociator, tracksTag, "RecoToTruthFixed"));
  trackAdaptiveToken_ = consumes<CaloMap>(mapTag(trackAssociator, tracksTag, "RecoToTruth" + adaptiveName_));

  for (auto const& pset : cfg.getParameter<std::vector<edm::ParameterSet>>("caloDomains")) {
    CaloDomainConfig domain;
    const auto name = pset.getParameter<std::string>("detector");
    if (name == "Ecal")
      domain.detector = Detector::Ecal;
    else if (name == "Hcal")
      domain.detector = Detector::Hcal;
    else if (name == "Hgcal")
      domain.detector = Detector::Hgcal;
    else
      throw cms::Exception("Configuration") << "caloDomains.detector must be Ecal, Hcal or Hgcal, got " << name;
    const auto collection = pset.getParameter<edm::InputTag>("collection");
    const auto associator = pset.getParameter<std::string>("associator");
    domain.isTrackster = pset.getParameter<bool>("tracksters");
    simEnergyScale_[static_cast<std::size_t>(domain.detector)] = pset.getParameter<double>("simEnergyScale");
    if (domain.isTrackster)
      domain.tracksters = consumes<std::vector<ticl::Trackster>>(collection);
    else
      domain.collection = consumes<edm::View<reco::CaloCluster>>(collection);
    domain.fixed = consumes<CaloMap>(mapTag(associator, collection, "RecoToTruthFixed"));
    domain.adaptive = consumes<CaloMap>(mapTag(associator, collection, "RecoToTruth" + adaptiveName_));
    domain.truthToReco = consumes<CaloMap>(mapTag(associator, collection, "TruthToReco"));
    caloDomains_.push_back(std::move(domain));
  }
  if (auto const tag = cfg.getParameter<edm::InputTag>("ticlCandidates"); !tag.label().empty()) {
    ticlCandidatesToken_ = consumes<std::vector<TICLCandidate>>(tag);
  }

  for (auto const& wp : workingPoints_)
    produces<CaloMap>(key_ + "RecoToTruth" + wp);
  produces<CaloMap>(key_ + "TruthToReco");
  produces<truth::PFCandidateTruthRecordCollection>(key_ + "TruthRecords");
  produces<truth::PFCandidateRecordCollection>(key_ + "CandidateRecords");
}

truth::CaloRegion PFCandidateTruthAssociator::regionOf(double absEta) const {
  if (absEta < barrelEtaMax_)
    return truth::CaloRegion::Barrel;
  if (absEta < endcapEtaMin_)
    return truth::CaloRegion::Transition;
  if (absEta < endcapEtaMax_)
    return truth::CaloRegion::Endcap;
  return truth::CaloRegion::Forward;
}

// Barrel candidates: walk the block elements. Endcap candidates (no blocks): the
// TICLCandidate recovered by index, whose tracksters are the constituents.
Constituents PFCandidateTruthAssociator::constituentsOf(reco::PFCandidate const& candidate,
                                                        edm::ProductID const& trackProductId,
                                                        std::vector<CaloDomain> const& domains,
                                                        std::vector<TICLCandidate> const* ticlCandidates,
                                                        std::vector<int32_t> const& ticlIndexOf,
                                                        std::size_t candidateIndex) const {
  Constituents out;
  if (candidate.trackRef().isNonnull() && candidate.trackRef().id() == trackProductId)
    out.track = static_cast<int32_t>(candidate.trackRef().key());

  auto slotFor = [&domains](edm::ProductID const& id) -> int {
    for (std::size_t s = 0; s < domains.size(); ++s)
      if (domains[s].productId == id)
        return static_cast<int>(s);
    return -1;
  };

  auto const& elements = candidate.elementsInBlocks();
  if (!elements.empty()) {
    for (auto const& [blockRef, index] : elements) {
      if (blockRef.isNull() || !blockRef.isAvailable())
        continue;
      auto const& element = blockRef->elements()[index];
      switch (element.type()) {
        case reco::PFBlockElement::ECAL:
        case reco::PFBlockElement::HCAL:
        case reco::PFBlockElement::HGCAL: {
          auto const& ref = element.clusterRef();
          if (ref.isNull())
            break;
          const int slot = slotFor(ref.id());
          if (slot >= 0)
            out.calo.emplace_back(static_cast<std::size_t>(slot), static_cast<uint32_t>(ref.key()));
          break;
        }
        case reco::PFBlockElement::TRACK:
          if (out.track < 0 && element.trackRef().isNonnull() && element.trackRef().id() == trackProductId)
            out.track = static_cast<int32_t>(element.trackRef().key());
          break;
        default:
          break;
      }
    }
    return out;
  }

  // No blocks: a pfTICL candidate. Its TICLCandidate is index-parallel in the pfTICL
  // collection, and the producer verified that mapping for this event.
  if (ticlCandidates == nullptr || candidateIndex >= ticlIndexOf.size() || ticlIndexOf[candidateIndex] < 0)
    return out;
  auto const& ticl = (*ticlCandidates)[ticlIndexOf[candidateIndex]];
  out.endcapMapped = true;
  for (auto const& ptr : ticl.tracksters()) {
    if (ptr.isNull())
      continue;
    const int slot = slotFor(ptr.id());
    if (slot >= 0)
      out.calo.emplace_back(static_cast<std::size_t>(slot), static_cast<uint32_t>(ptr.key()));
  }
  return out;
}

void PFCandidateTruthAssociator::produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const {
  auto const& graph = event.get(graphToken_);
  truth::SubgraphHitView hits(event.get(hitIndexToken_));
  auto const& targets = event.get(targetsToken_);
  auto const& candidates = event.get(candidatesToken_);
  const edm::Handle<std::vector<reco::Track>> tracksHandle = event.getHandle(tracksToken_);
  // The Fixed track map is consumed so the dependency is declared; the resolution reads
  // the adaptive rows only.
  event.get(trackFixedToken_);
  auto const& trackAdaptive = event.get(trackAdaptiveToken_).getMap();
  const uint32_t nParticles = graph.nParticles();
  const std::size_t nCandidates = candidates.size();

  // ---- constituent domains ------------------------------------------------------------
  std::vector<CaloDomain> domains;
  for (auto const& config : caloDomains_) {
    CaloDomain domain;
    domain.detector = config.detector;
    if (config.isTrackster) {
      const auto handle = event.getHandle(config.tracksters);
      if (!handle.isValid())
        continue;
      domain.productId = handle.id();
      domain.energies.reserve(handle->size());
      for (auto const& t : *handle)
        domain.energies.push_back(t.raw_energy());
    } else {
      const auto handle = event.getHandle(config.collection);
      if (!handle.isValid())
        continue;
      domain.productId = handle.id();
      domain.energies.reserve(handle->size());
      for (auto const& c : *handle)
        domain.energies.push_back(static_cast<float>(c.energy()));
    }
    const auto fixed = event.getHandle(config.fixed);
    const auto adaptive = event.getHandle(config.adaptive);
    const auto truthToReco = event.getHandle(config.truthToReco);
    if (fixed.isValid())
      domain.fixed = &fixed->getMap();
    if (adaptive.isValid())
      domain.adaptive = &adaptive->getMap();
    if (truthToReco.isValid())
      domain.truthToReco = &truthToReco->getMap();
    if (domain.valid())
      domains.push_back(std::move(domain));
  }

  // ---- endcap mapping: the k-th block-less candidate is pfTICL candidate k ------------
  std::vector<TICLCandidate> const* ticlCandidates = nullptr;
  std::vector<int32_t> ticlIndexOf(nCandidates, -1);
  if (!ticlCandidatesToken_.isUninitialized()) {
    const auto handle = event.getHandle(ticlCandidatesToken_);
    if (handle.isValid()) {
      std::vector<std::size_t> blockless;
      for (std::size_t i = 0; i < nCandidates; ++i)
        if (candidates[i].elementsInBlocks().empty())
          blockless.push_back(i);
      bool consistent = blockless.size() == handle->size();
      // Every charged pair must agree on the track, or the order is not what we assume.
      for (std::size_t k = 0; consistent && k < blockless.size(); ++k) {
        auto const& c = candidates[blockless[k]];
        auto const& t = (*handle)[k].trackPtr();
        const bool candidateCharged = c.trackRef().isNonnull();
        const bool ticlCharged = t.isNonnull();
        if (candidateCharged != ticlCharged || (candidateCharged && c.trackRef().key() != t.key()))
          consistent = false;
      }
      if (consistent) {
        ticlCandidates = &(*handle);
        for (std::size_t k = 0; k < blockless.size(); ++k)
          ticlIndexOf[blockless[k]] = static_cast<int32_t>(k);
      } else {
        std::call_once(endcapWarned_, [&] {
          edm::LogWarning("PFCandidateTruthAssociator")
              << "block-less candidates of '" << key_ << "' (" << blockless.size()
              << ") do not map index-parallel onto the TICLCandidates (" << handle->size()
              << "); endcap candidates get no calorimeter constituents this event";
        });
      }
    }
  }

  // ---- per-particle helpers -------------------------------------------------------------
  std::vector<bool> isTarget(nParticles, false);
  for (const unsigned int id : targets)
    if (id < nParticles)
      isTarget[id] = true;

  // Ancestor sets are needed for "is b a descendant of P" and for the merge test; cache
  // per particle on demand.
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> ancestorCache;
  auto ancestorsOf = [&](uint32_t id) -> std::unordered_set<uint32_t> const& {
    auto it = ancestorCache.find(id);
    if (it == ancestorCache.end()) {
      std::unordered_set<uint32_t> set;
      for (auto const& a : truth::Particle(&graph, id).ancestors())
        set.insert(a.id());
      it = ancestorCache.emplace(id, std::move(set)).first;
    }
    return it->second;
  };
  auto isDescendantOrSelf = [&](uint32_t b, uint32_t p) { return b == p || ancestorsOf(b).count(p) != 0u; };
  // The nearest validation-level particle at or above b, kNoBranch if none.
  auto targetAbove = [&](uint32_t b) -> uint32_t {
    if (b < nParticles && isTarget[b])
      return b;
    uint32_t best = kNoBranch;
    std::size_t bestDepth = std::numeric_limits<std::size_t>::max();
    for (const uint32_t a : ancestorsOf(b)) {
      if (isTarget[a]) {
        const std::size_t depth = ancestorsOf(a).size();
        if (best == kNoBranch || depth > bestDepth || (depth == bestDepth && a < best)) {
          // the deepest target ancestor is the nearest one
          best = a;
          bestDepth = depth;
        }
      }
    }
    return best;
  };
  auto hasTargetBelow = [&](uint32_t b) {
    for (auto const& d : truth::Particle(&graph, b).descendants())
      if (isTarget[d.id()])
        return true;
    return false;
  };

  // Sim energy per detector of a particle's subgraph.
  auto simEnergies = [&](uint32_t id) {
    std::array<double, kDetectors> e{};
    for (auto const& hit : hits.subgraphHits(truth::HitChannel::Calo, id))
      e[static_cast<std::size_t>(detectorOf(hit.detId))] += hit.energy;
    return e;
  };

  // Coverage of branch b by a set of calorimeter constituents, per detector: the sum of
  // the TruthToReco sharedEnergyFraction of those constituents in b's row.
  auto coverageOf = [&](uint32_t b, std::vector<std::pair<std::size_t, uint32_t>> const& calo) {
    std::array<double, kDetectors> cov{};
    for (auto const& [slot, index] : calo) {
      auto const& rows = *domains[slot].truthToReco;
      if (b >= rows.size())
        continue;
      for (auto const& e : rows[b])
        if (e.index() == index)
          cov[static_cast<std::size_t>(domains[slot].detector)] += e.sharedEnergy();
    }
    return cov;
  };

  // ---- track resolution -------------------------------------------------------------------
  // Best track per particle (the one sharing the most), from the adaptive rows.
  auto trackPurity = [&](std::size_t t) -> std::pair<uint32_t, double> {
    if (t >= trackAdaptive.size() || trackAdaptive[t].empty())
      return {kNoBranch, 0.};
    auto const& best = trackAdaptive[t][0];
    return {best.index(), 1. - static_cast<double>(best.score())};
  };
  std::vector<int32_t> trackOfParticle(nParticles, -1);
  std::vector<float> trackSharedOfParticle(nParticles, 0.f);
  for (std::size_t t = 0; t < trackAdaptive.size(); ++t) {
    const auto [p, purity] = trackPurity(t);
    if (p == kNoBranch || p >= nParticles || purity < minTrackPurity_)
      continue;
    const float shared = trackAdaptive[t][0].sharedEnergy();
    if (trackOfParticle[p] < 0 || shared > trackSharedOfParticle[p]) {
      trackOfParticle[p] = static_cast<int32_t>(t);
      trackSharedOfParticle[p] = shared;
    }
  }

  // ---- reco-driven: constituents, composition, level -----------------------------------------
  const edm::ProductID trackProductId = tracksHandle.isValid() ? tracksHandle.id() : edm::ProductID();
  std::vector<Constituents> constituents(nCandidates);
  std::vector<int32_t> candidateOfTrack(tracksHandle.isValid() ? tracksHandle->size() : 0, -1);
  // candidate -> (branch, attributed energy), the full composition
  std::vector<std::vector<std::pair<uint32_t, float>>> composition(nCandidates);
  auto records = std::make_unique<truth::PFCandidateRecordCollection>();
  records->reserve(nCandidates);

  for (std::size_t i = 0; i < nCandidates; ++i) {
    auto const& candidate = candidates[i];
    constituents[i] = constituentsOf(candidate, trackProductId, domains, ticlCandidates, ticlIndexOf, i);
    auto const& cons = constituents[i];
    if (cons.track >= 0 && static_cast<std::size_t>(cons.track) < candidateOfTrack.size())
      candidateOfTrack[cons.track] = static_cast<int32_t>(i);

    truth::PFCandidateRecord rec;
    rec.candidateIndex = static_cast<uint32_t>(i);
    rec.pfType = static_cast<int8_t>(candidate.particleId());
    rec.region = static_cast<uint8_t>(regionOf(std::abs(candidate.eta())));
    rec.energy = static_cast<float>(candidate.energy());
    rec.pt = static_cast<float>(candidate.pt());
    rec.eta = static_cast<float>(candidate.eta());
    rec.phi = static_cast<float>(candidate.phi());
    rec.hasTrack = cons.track >= 0;
    rec.endcapMapped = cons.endcapMapped;
    rec.nTracks = cons.track >= 0 ? 1 : 0;
    for (auto const& [slot, index] : cons.calo) {
      switch (domains[slot].detector) {
        case Detector::Ecal:
          ++rec.nEcal;
          break;
        case Detector::Hcal:
          ++rec.nHcal;
          break;
        default:
          ++rec.nTracksters;
      }
    }

    // Calorimeter composition: each constituent's energy goes to its adaptive best
    // branch; the whole-composition map gets every branch of its Fixed row weighted by
    // the constituent's purity for it.
    double caloEnergy = 0.;
    std::map<uint32_t, double> energyByBranch;  // adaptive attribution
    std::map<uint32_t, double> fixedByBranch;
    double unmatchedEnergy = 0.;
    for (auto const& [slot, index] : cons.calo) {
      auto const& domain = domains[slot];
      const double e = index < domain.energies.size() ? domain.energies[index] : 0.;
      caloEnergy += e;
      if (index < domain.adaptive->size() && !(*domain.adaptive)[index].empty()) {
        energyByBranch[(*domain.adaptive)[index][0].index()] += e;
      } else {
        unmatchedEnergy += e;
      }
      if (index < domain.fixed->size()) {
        for (auto const& el : (*domain.fixed)[index])
          fixedByBranch[el.index()] += e * std::max(0., 1. - static_cast<double>(el.score()));
      }
    }

    // Anchor: the track's particle.
    uint32_t anchor = kNoBranch;
    if (cons.track >= 0) {
      const auto [p, purity] = trackPurity(static_cast<std::size_t>(cons.track));
      if (p != kNoBranch && p < nParticles && purity >= minTrackPurity_) {
        anchor = p;
        rec.trackMatched = true;
      }
    }

    // Level.
    uint32_t level = kNoBranch;
    double own = 0., foreign = 0.;
    auto splitByRelation = [&](uint32_t base) {
      own = foreign = 0.;
      uint32_t dominantOther = kNoBranch;
      double dominantOtherEnergy = 0.;
      for (auto const& [b, e] : energyByBranch) {
        if (isDescendantOrSelf(b, base)) {
          own += e;
        } else {
          foreign += e;
          if (e > dominantOtherEnergy) {
            dominantOtherEnergy = e;
            dominantOther = b;
          }
        }
      }
      foreign += unmatchedEnergy;
      return dominantOther;
    };
    // Merge test: lift base to its common ancestor with other if the candidate covers
    // that ancestor well in some detector.
    auto tryMerge = [&](uint32_t base, uint32_t other) -> uint32_t {
      if (other == kNoBranch)
        return base;
      const auto lca = graph.lowestCommonAncestor({truth::Particle(&graph, base), truth::Particle(&graph, other)});
      if (!lca.has_value() || lca->id() >= nParticles || graph.particles()[lca->id()].isSynthetic())
        return base;
      const auto cov = coverageOf(lca->id(), cons.calo);
      const double best = *std::max_element(cov.begin(), cov.end());
      return best >= minMergedCoverage_ ? lca->id() : base;
    };

    if (anchor != kNoBranch) {
      level = anchor;
      const uint32_t other = splitByRelation(level);
      if (caloEnergy > 0. && foreign / caloEnergy > maxForeignFraction_) {
        const uint32_t merged = tryMerge(level, other);
        if (merged != level) {
          level = merged;
          splitByRelation(level);
        }
      }
    } else if (!energyByBranch.empty()) {
      // Calorimeter only: the dominant branch, then the merge test against the runner-up.
      uint32_t dominant = kNoBranch, runnerUp = kNoBranch;
      double e1 = 0., e2 = 0.;
      for (auto const& [b, e] : energyByBranch) {
        if (e > e1) {
          runnerUp = dominant;
          e2 = e1;
          dominant = b;
          e1 = e;
        } else if (e > e2) {
          runnerUp = b;
          e2 = e;
        }
      }
      level = dominant;
      splitByRelation(level);
      if (caloEnergy > 0. && foreign / caloEnergy > maxForeignFraction_) {
        const uint32_t merged = tryMerge(level, runnerUp);
        if (merged != level) {
          level = merged;
          splitByRelation(level);
        }
      }
    }

    rec.anchorId = anchor;
    rec.branchId = level;
    if (caloEnergy > 0.) {
      rec.ownFraction = static_cast<float>(own / caloEnergy);
      rec.foreignFraction = static_cast<float>(foreign / caloEnergy);
    } else if (level != kNoBranch) {
      rec.ownFraction = 1.f;  // track-only
    }

    // Composition for the maps: the level with the candidate's own energy, plus every
    // other branch of the Fixed composition.
    if (level != kNoBranch) {
      const double levelEnergy = caloEnergy > 0. ? own : candidate.energy();
      composition[i].emplace_back(level, static_cast<float>(levelEnergy));
      for (auto const& [b, e] : fixedByBranch)
        if (b != level)
          composition[i].emplace_back(b, static_cast<float>(e));
    }
    records->push_back(rec);
  }

  // ---- reco-driven classes -----------------------------------------------------------------
  // The validation-level particle each candidate is judged against; a target claimed by
  // several candidates keeps the highest-energy one as Matched and the others as Split.
  std::unordered_map<uint32_t, std::size_t> claimant;
  std::vector<uint32_t> targetOf(nCandidates, kNoBranch);
  for (std::size_t i = 0; i < nCandidates; ++i) {
    auto& rec = (*records)[i];
    if (rec.branchId == kNoBranch) {
      const bool noDomainConstituent = constituents[i].track < 0 && constituents[i].calo.empty();
      rec.candidateClass = static_cast<uint8_t>(noDomainConstituent ? PFCandidateClass::Unevaluated
                                                                     : PFCandidateClass::Unmatched);
      continue;
    }
    const uint32_t target = targetAbove(rec.branchId);
    if (target == kNoBranch) {
      rec.candidateClass = static_cast<uint8_t>(hasTargetBelow(rec.branchId) ? PFCandidateClass::Merged
                                                                              : PFCandidateClass::Other);
      continue;
    }
    rec.targetId = target;
    targetOf[i] = target;
    auto it = claimant.find(target);
    if (it == claimant.end() || candidates[i].energy() > candidates[it->second].energy())
      claimant[target] = i;
  }
  for (std::size_t i = 0; i < nCandidates; ++i) {
    auto& rec = (*records)[i];
    if (targetOf[i] == kNoBranch)
      continue;
    rec.candidateClass = static_cast<uint8_t>(claimant[targetOf[i]] == i ? PFCandidateClass::Matched
                                                                          : PFCandidateClass::Split);
  }

  // ---- maps ----------------------------------------------------------------------------------
  auto truthToReco = std::make_unique<CaloMap>(nParticles);
  for (auto const& wp : workingPoints_) {
    auto map = std::make_unique<CaloMap>(nCandidates);
    const bool full = wp == "Fixed";
    for (std::size_t i = 0; i < nCandidates; ++i) {
      double total = 0.;
      for (auto const& [b, e] : composition[i])
        total += e;
      for (std::size_t k = 0; k < composition[i].size(); ++k) {
        auto const& [b, e] = composition[i][k];
        const float score = total > 0. ? static_cast<float>(1. - e / total) : 1.f;
        map->insert(static_cast<unsigned int>(i), b, e, score);
        if (!full)
          break;  // the level only
      }
    }
    map->sort(byAscendingScore);
    event.put(std::move(map), key_ + "RecoToTruth" + wp);
  }
  for (std::size_t i = 0; i < nCandidates; ++i) {
    if (composition[i].empty())
      continue;
    auto const& [level, e] = composition[i][0];
    const double truthEnergy = graph.particles()[level].momentum.energy();
    const float fraction = truthEnergy > 0. ? static_cast<float>(std::min(1., e / truthEnergy)) : 0.f;
    truthToReco->insert(level, static_cast<unsigned int>(i), fraction, 1.f - fraction);
  }
  truthToReco->sort(byAscendingScore);
  event.put(std::move(truthToReco), key_ + "TruthToReco");

  // ---- truth-driven ladder ---------------------------------------------------------------------
  auto ladder = std::make_unique<truth::PFCandidateTruthRecordCollection>();
  ladder->reserve(targets.size());
  for (const unsigned int p : targets) {
    if (p >= nParticles)
      continue;
    auto const& data = graph.particles()[p];
    truth::PFCandidateTruthRecord r;
    r.particleId = p;
    r.pdgId = data.pdgId;
    r.energy = static_cast<float>(data.momentum.energy());
    r.pt = static_cast<float>(data.momentum.pt());
    r.eta = static_cast<float>(data.momentum.eta());
    r.phi = static_cast<float>(data.momentum.phi());
    const auto entry = truth::Particle(&graph, p).checkpoint(0);
    r.caloEta = entry.has_value() ? static_cast<float>(entry->position.eta()) : r.eta;
    const auto region = regionOf(std::abs(r.caloEta));
    r.region = static_cast<uint8_t>(region);

    const auto sim = simEnergies(p);
    r.simEnergyEcal = static_cast<float>(sim[0]);
    r.simEnergyHcal = static_cast<float>(sim[1]);
    r.simEnergyHgcal = static_cast<float>(sim[2]);
    // Expected detectors: a share of the particle's calorimeter sim energy above the
    // floor, and compatible with the region. Sampling energies differ between detectors,
    // so the share is taken within the detectors of the region only.
    std::array<bool, kDetectors> expected{};
    {
      std::array<bool, kDetectors> inRegion{};
      inRegion[0] = inRegion[1] = region != truth::CaloRegion::Endcap && region != truth::CaloRegion::Forward;
      inRegion[2] = region != truth::CaloRegion::Barrel;
      std::array<double, kDetectors> scaled{};
      double total = 0.;
      for (std::size_t d = 0; d < kDetectors; ++d) {
        scaled[d] = sim[d] * simEnergyScale_[d];
        if (inRegion[d])
          total += scaled[d];
      }
      for (std::size_t d = 0; d < kDetectors; ++d)
        expected[d] = inRegion[d] && total > 0. && scaled[d] / total >= minExpectedDetectorFraction_;
    }
    r.set(PFRung::EcalExpected, expected[0]);
    r.set(PFRung::HcalExpected, expected[1]);
    r.set(PFRung::HgcalExpected, expected[2]);

    // Track rungs.
    const bool charged = isCharged(data.pdgId);
    r.set(PFRung::TrackExpected, charged);
    if (charged && trackOfParticle[p] >= 0) {
      r.set(PFRung::TrackFound);
      r.trackIndex = trackOfParticle[p];
      if (static_cast<std::size_t>(r.trackIndex) < candidateOfTrack.size() && candidateOfTrack[r.trackIndex] >= 0) {
        r.set(PFRung::TrackInCandidate);
        r.candidateIndex = candidateOfTrack[r.trackIndex];
      }
    }

    // Collected: per detector, the summed sharedEnergyFraction over all clusters in the
    // particle's TruthToReco row. Also remember which candidate holds most of it.
    std::array<double, kDetectors> collected{};
    std::array<std::map<int32_t, double>, kDetectors> byCandidate;  // candidate -> collected share
    // constituent (slot, index) -> candidate holding it
    std::unordered_map<uint64_t, int32_t> holder;
    for (std::size_t i = 0; i < nCandidates; ++i)
      for (auto const& [slot, index] : constituents[i].calo)
        holder[(static_cast<uint64_t>(slot) << 32) | index] = static_cast<int32_t>(i);
    for (std::size_t slot = 0; slot < domains.size(); ++slot) {
      auto const& rows = *domains[slot].truthToReco;
      if (p >= rows.size())
        continue;
      const auto d = static_cast<std::size_t>(domains[slot].detector);
      for (auto const& e : rows[p]) {
        collected[d] += e.sharedEnergy();
        const auto it = holder.find((static_cast<uint64_t>(slot) << 32) | e.index());
        if (it != holder.end())
          byCandidate[d][it->second] += e.sharedEnergy();
      }
    }
    r.collectedEcal = static_cast<float>(std::min(1., collected[0]));
    r.collectedHcal = static_cast<float>(std::min(1., collected[1]));
    r.collectedHgcal = static_cast<float>(std::min(1., collected[2]));
    r.set(PFRung::EcalCollected, expected[0] && collected[0] >= minCollectedFraction_);
    r.set(PFRung::HcalCollected, expected[1] && collected[1] >= minCollectedFraction_);
    r.set(PFRung::HgcalCollected, expected[2] && collected[2] >= minCollectedFraction_);

    // Resolved candidate without a track: the one holding most collected energy, summed
    // over the expected detectors.
    if (r.candidateIndex == truth::PFCandidateTruthRecord::kNone) {
      std::map<int32_t, double> total;
      for (std::size_t d = 0; d < kDetectors; ++d)
        if (expected[d])
          for (auto const& [c, share] : byCandidate[d])
            total[c] += share;
      int32_t best = truth::PFCandidateTruthRecord::kNone;
      double bestShare = 0.;
      for (auto const& [c, share] : total)
        if (share > bestShare) {
          bestShare = share;
          best = c;
        }
      r.candidateIndex = best;
    }

    // Linked: the share of the collected energy that sits in the resolved candidate.
    const int32_t K = r.candidateIndex;
    auto linked = [&](std::size_t d) {
      if (K == truth::PFCandidateTruthRecord::kNone || collected[d] <= 0.)
        return 0.;
      const auto it = byCandidate[d].find(K);
      return it == byCandidate[d].end() ? 0. : std::min(1., it->second / collected[d]);
    };
    r.linkedEcal = static_cast<float>(linked(0));
    r.linkedHcal = static_cast<float>(linked(1));
    r.linkedHgcal = static_cast<float>(linked(2));
    r.set(PFRung::EcalLinked, r.has(PFRung::EcalCollected) && r.linkedEcal >= minLinkedFraction_);
    r.set(PFRung::HcalLinked, r.has(PFRung::HcalCollected) && r.linkedHcal >= minLinkedFraction_);
    r.set(PFRung::HgcalLinked, r.has(PFRung::HgcalCollected) && r.linkedHgcal >= minLinkedFraction_);

    // Calo-to-calo: the candidate holding most of the ECAL piece is the one holding most
    // of the HCAL piece.
    if (r.has(PFRung::EcalCollected) && r.has(PFRung::HcalCollected)) {
      auto holderOf = [&](std::size_t d) {
        int32_t best = -1;
        double bestShare = 0.;
        for (auto const& [c, share] : byCandidate[d])
          if (share > bestShare) {
            bestShare = share;
            best = c;
          }
        return best;
      };
      const int32_t e = holderOf(0), h = holderOf(1);
      r.set(PFRung::CaloSameCandidate, e >= 0 && e == h);
    }

    // Candidate rungs, judged on the resolved candidate.
    if (K != truth::PFCandidateTruthRecord::kNone) {
      auto const& crec = (*records)[K];
      r.candidateType = crec.pfType;
      r.foreignFraction = crec.foreignFraction;
      r.energyRatio = r.energy > 0.f ? crec.energy / r.energy : 0.f;
      if (crec.branchId != kNoBranch && isDescendantOrSelf(p, crec.branchId)) {
        r.set(PFRung::CandidateFound);
        r.set(PFRung::Merged, crec.branchId != p);
      }
      r.set(PFRung::Clean, r.has(PFRung::CandidateFound) && crec.foreignFraction <= maxForeignFraction_);
      const int expectedType = expectedPFType(data.pdgId);
      r.expectedType = static_cast<int8_t>(expectedType);
      if (expectedType >= 0) {
        r.set(PFRung::PdgEvaluated);
        r.set(PFRung::PdgCorrect, r.has(PFRung::CandidateFound) && crec.pfType == expectedType);
      }
    } else {
      const int expectedType = expectedPFType(data.pdgId);
      r.expectedType = static_cast<int8_t>(expectedType);
      r.set(PFRung::PdgEvaluated, expectedType >= 0);
    }
    ladder->push_back(r);
  }

  event.put(std::move(ladder), key_ + "TruthRecords");
  event.put(std::move(records), key_ + "CandidateRecords");
}

void PFCandidateTruthAssociator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<edm::InputTag>("targets", edm::InputTag("truthBranchTargets", "truthToRecoTargetsReconstructableFromSignal"))
      ->setComment("the validation level: one record per particle of this list");
  desc.add<edm::InputTag>("candidates", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("trackAssociator", "allTrackToTruthBranchAssociators");
  desc.add<edm::InputTag>("ticlCandidates", edm::InputTag("ticlCandidate"))
      ->setComment("TICLCandidates of the pfTICL candidates, recovered index-parallel; empty label disables the endcap");
  edm::ParameterSetDescription calo;
  calo.add<std::string>("detector")->setComment("Ecal, Hcal or Hgcal");
  calo.add<edm::InputTag>("collection");
  calo.add<std::string>("associator");
  calo.add<bool>("tracksters", false);
  calo.add<double>("simEnergyScale", 1.)
      ->setComment(
          "reco energy per unit of sim-hit energy in this detector, for the expected-detector rule. Measured on "
          "no-PU ttbar D110 from clusters fully covered by one branch: ECAL 1.05, HCAL 105, HGCAL taken as 1");
  desc.addVPSet("caloDomains",
                calo,
                {[] {
                   edm::ParameterSet p;
                   p.addParameter<std::string>("detector", "Ecal");
                   p.addParameter<edm::InputTag>("collection", edm::InputTag("particleFlowClusterECAL"));
                   p.addParameter<std::string>("associator", "truthBranchPFClusterEcalAssociators");
                   p.addParameter<bool>("tracksters", false);
                   p.addParameter<double>("simEnergyScale", 1.05);
                   return p;
                 }(),
                 [] {
                   edm::ParameterSet p;
                   p.addParameter<std::string>("detector", "Hcal");
                   p.addParameter<edm::InputTag>("collection", edm::InputTag("particleFlowClusterHCAL"));
                   p.addParameter<std::string>("associator", "truthBranchPFClusterHcalAssociators");
                   p.addParameter<bool>("tracksters", false);
                   p.addParameter<double>("simEnergyScale", 105.);
                   return p;
                 }(),
                 [] {
                   edm::ParameterSet p;
                   p.addParameter<std::string>("detector", "Hgcal");
                   p.addParameter<edm::InputTag>("collection", edm::InputTag("ticlCandidate"));
                   p.addParameter<std::string>("associator", "truthBranchTracksterAssociators");
                   p.addParameter<bool>("tracksters", true);
                   p.addParameter<double>("simEnergyScale", 1.);
                   return p;
                 }()});
  desc.add<std::vector<std::string>>("workingPointNames", {"Fixed", "AdaptiveTight", "AdaptiveNominal"})
      ->setComment("RecoToTruth products: Fixed carries the whole composition, every other name the resolved level");
  desc.add<std::string>("constituentWorkingPoint", "AdaptiveNominal")
      ->setComment("the constituent maps' working point that attributes each constituent to one branch");
  desc.add<double>("minTrackPurity", 0.75);
  desc.add<double>("minCollectedFraction", 0.5);
  desc.add<double>("minLinkedFraction", 0.5);
  desc.add<double>("maxForeignFraction", 0.25);
  desc.add<double>("minExpectedDetectorFraction", 0.1);
  desc.add<double>("minMergedCoverage", 0.5)
      ->setComment("a common ancestor becomes the level only if the candidate covers this share of it in some detector");
  desc.add<double>("barrelEtaMax", 1.48);
  desc.add<double>("endcapEtaMin", 1.6);
  desc.add<double>("endcapEtaMax", 3.0);
  descriptions.add("pfCandidateTruthAssociator", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCandidateTruthAssociator);
