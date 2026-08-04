#ifndef SimDataFormats_TrackingAnalysis_SimNtuplet_h
#define SimDataFormats_TrackingAnalysis_SimNtuplet_h

#include "SimPixelTrack.h"

/** @brief SimPixelTrack::Ntuplet or short SimNtuplet is a true RecHit N-tuplet 
 *         of arbitrary length of a simulated particle/TrackingParticle.
 *
 * SimNtuplets are defined to be chains of SimDoublets that pairwise share the middle RecHit.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date April 2026
 */

class SimPixelTrack::Ntuplet {
public:
  // flags indicating qualities of Ntuplet (depending on its constituents)
  // The order is chosen in such a way that a smaller status value means that the Ntuplet get farther
  // in the reconstruction chain. Hence, a value of 0 corresponds to the Ntuplet surviving reconstruction.
  enum class StatusBit : status_type {
    isTooShort = 1,
    hasMissingLayerPair = 1 << 1,
    hasUndefDoubletCuts = 1 << 2,
    hasKilledDoublets = 1 << 3,
    hasUndefTripletCuts = 1 << 4,
    hasKilledTriplets = 1 << 5,
    hasKilledQuadruplets = 1 << 6,
    invalidStart = 1 << 7
  };

  // default constructor
  Ntuplet() = default;

  // constructor
  Ntuplet(size_t const numDoublets,
          status_type const status,
          layer_type const firstLayerId,
          layer_type const secondLayerId,
          layer_type const lastLayerId,
          int_type const numSkippedLayers)
      : numDoublets_(numDoublets),
        status_(status),
        firstLayerId_(firstLayerId),
        secondLayerId_(secondLayerId),
        lastLayerId_(lastLayerId),
        numSkippedLayers_(numSkippedLayers) {};

  // accessing the different members
  size_t numDoublets() const { return numDoublets_; }
  size_t numRecHits() const { return (numDoublets_ + 1); }
  layer_type firstLayerId() const { return firstLayerId_; }
  layer_type secondLayerId() const { return secondLayerId_; }
  layer_type lastLayerId() const { return lastLayerId_; }
  int_type numSkippedLayers() const { return numSkippedLayers_; }

  // method to update an external status
  static status_type updateStatus(status_type status,
                                  bool const hasUndefDoubletCuts,
                                  bool const hasMissingLayerPair,
                                  bool const hasKilledDoublets,
                                  bool const hasUndefTripletCuts,
                                  bool const hasKilledTriplets,
                                  bool const hasKilledQuadruplets,
                                  bool const isTooShort = false,
                                  bool const invalidStart = false) {
    return status | (status_type(hasUndefDoubletCuts) * status_type(StatusBit::hasUndefDoubletCuts) +
                     status_type(hasMissingLayerPair) * status_type(StatusBit::hasMissingLayerPair) +
                     status_type(hasKilledDoublets) * status_type(StatusBit::hasKilledDoublets) +
                     status_type(hasUndefTripletCuts) * status_type(StatusBit::hasUndefTripletCuts) +
                     status_type(hasKilledTriplets) * status_type(StatusBit::hasKilledTriplets) +
                     status_type(hasKilledQuadruplets) * status_type(StatusBit::hasKilledQuadruplets) +
                     status_type(isTooShort) * status_type(StatusBit::isTooShort) +
                     status_type(invalidStart) * status_type(StatusBit::invalidStart));
  }

  // methods to set status to alive, undef or killed
  void setUndefDoubletCuts() { status_ |= status_type(StatusBit::hasUndefDoubletCuts); }
  void setUndefTripletCuts() { status_ |= status_type(StatusBit::hasUndefTripletCuts); }
  void setMissingLayerPair() { status_ |= status_type(StatusBit::hasMissingLayerPair); }
  void setKilledDoublets() { status_ |= status_type(StatusBit::hasKilledDoublets); }
  void setKilledTriplet() { status_ |= status_type(StatusBit::hasKilledTriplets); }
  void setKilledQuadruplet() { status_ |= status_type(StatusBit::hasKilledQuadruplets); }
  void setTooShort() { status_ |= status_type(StatusBit::isTooShort); }
  void setInvalidStart() { status_ |= status_type(StatusBit::invalidStart); }

  // methods to check if status is undef, alive or killed
  bool hasUndefDoubletCuts() const { return status_ & status_type(StatusBit::hasUndefDoubletCuts); }
  bool hasUndefTripletCuts() const { return status_ & status_type(StatusBit::hasUndefTripletCuts); }
  bool hasUndef() const { return hasUndefDoubletCuts() || hasUndefTripletCuts(); }
  bool hasMissingLayerPair() const { return status_ & status_type(StatusBit::hasMissingLayerPair); }
  bool hasKilledDoublets() const { return status_ & status_type(StatusBit::hasKilledDoublets); }
  bool hasKilledTriplets() const { return status_ & status_type(StatusBit::hasKilledTriplets); }
  bool hasKilledQuadruplets() const { return status_ & status_type(StatusBit::hasKilledQuadruplets); }
  bool isKilled() const {
    return hasMissingLayerPair() || hasKilledDoublets() || hasKilledTriplets() || hasKilledQuadruplets();
  }
  bool isTooShort() const { return status_ & status_type(StatusBit::isTooShort); }
  bool invalidStart() const { return status_ & status_type(StatusBit::invalidStart); }
  bool isAlive() const { return !(status_); }  // if nothing is set (no undef and no kills) the tuplet is alive

  // method to get the leading digit of the status (first non-zero one),
  // e.g. status=00110100 -> failingRecoStep()=00000100
  // This represents the first step of the reco chain the given Ntuplet fails.
  // For an alive Ntuplet, return the max value 11111111.
  status_type failingRecoStep() const {
    if (isAlive())
      return 0b11111111;
    else
      return (status_ & ((~status_) + 1));
  }

  // method to compare the own status to a given reference and check which one gets farther in the reconstruction chain
  bool getsFartherInRecoChainThanReference(Ntuplet const& referenceNtuplet) const {
    return failingRecoStep() > referenceNtuplet.failingRecoStep();
  }
  // method to check if the own Ntuplet gets exactly as far in the reco chain as a given reference
  bool getsAsFarInRecoChainAsReference(Ntuplet const& referenceNtuplet) const {
    return failingRecoStep() == referenceNtuplet.failingRecoStep();
  }

private:
  size_t numDoublets_;       // number of doublets in the Ntuplet
  status_type status_;       // status flags of the Ntuplet (missing layer pairs, undefined cuts, killed doublets, etc.)
  layer_type firstLayerId_;  // index of the first layer of the Ntuplet
  layer_type secondLayerId_;   // index of the second layer of the Ntuplet
  layer_type lastLayerId_;     // index of the last layer of the Ntuplet
  int_type numSkippedLayers_;  // number of skipped layers over the full Ntuplet (sum of skips by doublets)
};

#endif
