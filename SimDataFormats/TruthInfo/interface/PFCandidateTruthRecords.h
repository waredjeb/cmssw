#ifndef SimDataFormats_TruthInfo_interface_PFCandidateTruthRecords_h
#define SimDataFormats_TruthInfo_interface_PFCandidateTruthRecords_h

// Persisted records of the particle-flow candidate truth validation.
//
// Two records, both filled by PFCandidateTruthAssociator from the constituent
// association maps (tracks, PFClusters, tracksters) and the truth graph:
//
//   PFCandidateTruthRecord   one per truth particle of the validation level (the
//                            reconstructableFromSignal antichain): the findings of the
//                            reconstruction ladder, stored UNCONDITIONALLY as flags so a
//                            consumer can draw them as a cut-flow or one rung at a time.
//   PFCandidateRecord        one per reco::PFCandidate: the truth branch it resolves to
//                            through the track-first rule, its class (matched, merged,
//                            split, fake) and its truth composition.
//
// The thresholds that decide each flag live on the producer; the records carry the
// underlying fractions too, so a threshold can be revisited offline.

#include <cstdint>
#include <limits>
#include <vector>

namespace truth {

  // Where a particle or candidate is judged, from the truth eta at the calorimeter
  // entrance (particles) or the candidate direction (candidates). The transition bin is
  // kept apart because both calorimeter systems can be expected there.
  enum class CaloRegion : uint8_t { Barrel = 0, Transition = 1, Endcap = 2, Forward = 3 };

  // One bit per rung of the reconstruction ladder. "Expected" bits say whether the rung
  // applies to this particle at all: a rung that is not expected is neither passed nor
  // failed. A consumer must gate every rung on its Expected bit.
  enum class PFRung : uint32_t {
    TrackExpected = 1u << 0,      // charged particle
    TrackFound = 1u << 1,         // a track resolves to this particle above the purity floor
    TrackInCandidate = 1u << 2,   // that track is the track of some candidate
    EcalExpected = 1u << 3,       // sim energy in ECAL above the expected-detector floor
    EcalCollected = 1u << 4,      // ECAL clusters cover the particle's ECAL deposit
    EcalLinked = 1u << 5,         // those clusters are in the particle's candidate
    HcalExpected = 1u << 6,
    HcalCollected = 1u << 7,
    HcalLinked = 1u << 8,
    HgcalExpected = 1u << 9,
    HgcalCollected = 1u << 10,
    HgcalLinked = 1u << 11,
    CaloSameCandidate = 1u << 12,  // the ECAL and HCAL pieces sit in one candidate
    CandidateFound = 1u << 13,     // a candidate resolves to this particle or an ancestor
    Merged = 1u << 14,             // that candidate resolves to an ancestor
    Clean = 1u << 15,              // foreign energy in the candidate below the floor
    PdgEvaluated = 1u << 16,       // the species has an expected PF type
    PdgCorrect = 1u << 17
  };

  [[nodiscard]] constexpr uint32_t rungBit(PFRung rung) { return static_cast<uint32_t>(rung); }

  struct PFCandidateTruthRecord {
    static constexpr int32_t kNone = -1;

    uint32_t particleId = 0;  // index in truth::Graph
    int32_t pdgId = 0;
    uint8_t region = static_cast<uint8_t>(CaloRegion::Barrel);
    uint32_t rungs = 0;  // PFRung bits

    // The candidate the particle resolves to: the one holding its track when a track was
    // found, else the one holding most of its collected calorimeter energy. kNone if none.
    int32_t candidateIndex = kNone;
    int32_t trackIndex = kNone;  // the track found for it, kNone if none
    int8_t candidateType = -1;   // reco::PFCandidate::ParticleType of that candidate
    int8_t expectedType = -1;    // the species' expected PF type, -1 if not evaluated

    // Truth kinematics: the generator or SimTrack momentum, and the eta at the
    // calorimeter entrance the region is decided from.
    float energy = 0.f;
    float pt = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    float caloEta = 0.f;

    // Sim energy of the particle's subgraph per detector, from the hit index, as stored
    // there: sampling energies, comparable within one detector only. The expected-detector
    // rule scales them per detector to a common energy scale before comparing.
    float simEnergyEcal = 0.f;
    float simEnergyHcal = 0.f;
    float simEnergyHgcal = 0.f;

    // Fraction of the particle's deposit in each detector held by any cluster
    // (collected), and the share of that collected energy inside the resolved
    // candidate (linked). Both in [0, 1].
    float collectedEcal = 0.f;
    float collectedHcal = 0.f;
    float collectedHgcal = 0.f;
    float linkedEcal = 0.f;
    float linkedHcal = 0.f;
    float linkedHgcal = 0.f;

    // Of the resolved candidate: the energy share carried by constituents that belong
    // neither to this particle nor to its descendants, and its energy over the truth
    // energy. The response is a distribution, never a gate.
    float foreignFraction = 0.f;
    float energyRatio = 0.f;

    [[nodiscard]] bool has(PFRung rung) const { return (rungs & rungBit(rung)) != 0u; }
    void set(PFRung rung, bool value = true) {
      if (value)
        rungs |= rungBit(rung);
      else
        rungs &= ~rungBit(rung);
    }
  };

  enum class PFCandidateClass : uint8_t {
    // No constituent resolves to any SELECTED branch. The constituent associators
    // match against the branch selector's candidates (a pt floor on stable particles),
    // so a candidate built on a sub-threshold particle lands here next to a genuine
    // fake; the two are told apart only by lowering the floor.
    Unmatched = 0,
    Matched = 1,     // resolves to a validation-level particle no other candidate claims
    Merged = 2,      // resolves to an ancestor of validation-level particles
    Split = 3,       // resolves to a particle another candidate already claims
    Other = 4,       // resolves to a branch outside the validation level (e.g. pileup)
    Unevaluated = 5  // no constituent in any configured domain (e.g. HF-only candidates)
  };

  struct PFCandidateRecord {
    static constexpr uint32_t kNoBranch = std::numeric_limits<uint32_t>::max();

    uint32_t candidateIndex = 0;
    int8_t pfType = -1;  // reco::PFCandidate::ParticleType
    uint8_t region = static_cast<uint8_t>(CaloRegion::Barrel);
    uint8_t candidateClass = static_cast<uint8_t>(PFCandidateClass::Unmatched);
    bool hasTrack = false;
    bool trackMatched = false;  // the track resolved to a branch above the purity floor
    bool endcapMapped = false;  // pfTICL candidate whose TICLCandidate was recovered

    // The branch the candidate resolves to after the track-first rule (the level), the
    // particle its track resolved to (the anchor, kNoBranch without a track), and the
    // validation-level particle the class is decided against.
    uint32_t branchId = kNoBranch;
    uint32_t anchorId = kNoBranch;
    uint32_t targetId = kNoBranch;

    float energy = 0.f;
    float pt = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    // Energy share of the calorimeter constituents that belong to the level branch or
    // its descendants, and the share that belongs to no related branch.
    float ownFraction = 0.f;
    float foreignFraction = 0.f;
    uint8_t nTracks = 0;
    uint8_t nEcal = 0;
    uint8_t nHcal = 0;
    uint8_t nTracksters = 0;
  };

  using PFCandidateTruthRecordCollection = std::vector<PFCandidateTruthRecord>;
  using PFCandidateRecordCollection = std::vector<PFCandidateRecord>;

}  // namespace truth

#endif
