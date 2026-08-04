#ifndef SimDataFormats_TrackingAnalysis_SimTriplet_h
#define SimDataFormats_TrackingAnalysis_SimTriplet_h

#include "SimPixelTrack.h"

/** @brief SimPixelTrack::Triplet or short SimTriplet is a true (shortest) RecHit triplet 
 *         of a simulated particle/TrackingParticle.
 *
 * SimTriplets are defined to be pairs of SimDoublets that share the middle RecHit.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date April 2026
 */

class SimPixelTrack::Triplet {
public:
  enum class Status : status_type { undef, alive, killedByCuts, killedByMissingLayerPair };
  
  Triplet(size_t innerDoubletIndex, size_t nInnerTriplets)
      : index_(innerDoubletIndex), status_(Status::undef), quadrupletIsKilled_(nInnerTriplets, false) {}

  size_t innerDoubletIndex() const { return index_; }

  // methods to set status to undef, alive or killed
  void setUndef() { status_ = Status::undef; }
  void setAlive() { status_ = Status::alive; }
  void setKilled() { status_ = Status::killedByCuts; }

  // methods to check if status is undef, alive or killed
  bool isUndef() const { return status_ == Status::undef; }
  bool isAlive() const { return status_ == Status::alive; }
  bool isKilled() const { return status_ == Status::killedByCuts; }

  void setCurvature(float_type const curvature) { curvature_ = curvature; }
  float_type curvature() const { return curvature_; }

  void setKilledQuadruplet(size_t i) { quadrupletIsKilled_.at(i) = true; }
  bool isKilledQuadruplet(size_t i) { return quadrupletIsKilled_.at(i); }
  std::vector<bool> const& quadruplets() const { return quadrupletIsKilled_; }

private:
  size_t index_;                            // index of the inner doublet of the triplet
  Status status_;                           // status of the triplet
  float_type curvature_{-99999};            // curvature of the triplet
  std::vector<bool> quadrupletIsKilled_{};  // status of the quadruplets with this triplet as the the outer triplet
};

#endif
