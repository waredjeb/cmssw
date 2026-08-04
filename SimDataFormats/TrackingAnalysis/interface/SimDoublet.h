#ifndef SimDataFormats_TrackingAnalysis_SimDoublet_h
#define SimDataFormats_TrackingAnalysis_SimDoublet_h

#include "SimPixelTrack.h"

/** @brief SimPixelTrack::Doublet or short SimDoublet is a true (shortest) RecHit doublet 
 *         of a simulated particle/TrackingParticle.
 *
 * SimPixelTracks hold references to all pixel RecHits of a simulated TrackingParticle.
 * Ones those RecHits are sorted according to their position relative to the particle vertex
 * by the method sortRecHits(), you can create the true doublets of RecHits that the 
 * TrackingParticle left in the detector. These SimPixelTrack::Doublet objects can be used  
 * to optimize the doublet creation in the reconstruction.
 *
 * The Doublets are generated as the RecHit pairs between two consecutively hit layers.
 * I.e., if a TrackingParticle produces
 *  - 1 hit (A) in 1st layer
 *  - 2 hits (B, C) in 3rd layer
 *  - 1 hit (D) in 4th layer
 * then, the true Doublets are:
 *  (A-B), (A-C), (B-D) and (C-D).
 * So, neither does it matter that the 2nd layer got "skipped" as there are no hits,
 * nor is the Doublet of (A-D) formed since there is a layer with hits in between.
 * Doublets are not created between hits within the same layer.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date April 2026
 */

class SimPixelTrack::Doublet {
public:
  // possible states of the doublet (could be set by an analyzer according to doublet cuts)
  enum class Status : status_type { undef, alive, killedByCuts, killedByMissingLayerPair };

  // constructors
  Doublet() = delete;
  Doublet(SimPixelTrack const&, size_t const, size_t const, const TrackerTopology*, std::vector<size_t> const&);

  // destructor
  ~Doublet();

  Doublet& operator=(const Doublet& other) {
    if (this != &other) {
      // Note: the reference to the SimPixelTrack will remain unchanged!
      innerHit_ = other.innerHit_;
      outerHit_ = other.outerHit_;
    }
    return *this;
  }

  // method to get the number of skipped layers
  int_type numSkippedLayers() const { return numSkippedLayers_; }

  // method to get the layer pair ID
  layer_type layerPairId() const { return layerPairId_; }

  // methods to get the inner/outer layerId
  layer_type innerLayerId() const { return simPixelTrack_.hits_.layerIds[innerHit_]; }
  layer_type outerLayerId() const { return simPixelTrack_.hits_.layerIds[outerHit_]; }

  // methods to get the cluster size of the inner/outer RecHit
  int_type innerClusterYSize() const { return simPixelTrack_.hits_.clusterYSizes[innerHit_]; }
  int_type outerClusterYSize() const { return simPixelTrack_.hits_.clusterYSizes[outerHit_]; }

  // methods to get the module ids of the inner/outer RecHit
  unsigned int innerModuleId() const { return simPixelTrack_.hits_.moduleIds[innerHit_]; }
  unsigned int outerModuleId() const { return simPixelTrack_.hits_.moduleIds[outerHit_]; }

  // methods to get the global position of the inner/outer RecHit
  const GlobalPoint& innerGlobalPos() const { return simPixelTrack_.hits_.globalPositions[innerHit_]; };
  const GlobalPoint& outerGlobalPos() const { return simPixelTrack_.hits_.globalPositions[outerHit_]; };

  // methods to set status to undef, alive or killed
  void setUndef() { status_ = Status::undef; }
  void setAlive() { status_ = Status::alive; }
  void setKilledByCuts() { status_ = Status::killedByCuts; }
  void setKilledByMissingLayerPair() { status_ = Status::killedByMissingLayerPair; }
  void setValidStart() { validStart_ = true; }

  // methods to check if status is undef, alive or killed
  bool isUndef() const { return status_ == Status::undef; }
  bool isAlive() const { return status_ == Status::alive; }
  bool isKilledByCuts() const { return status_ == Status::killedByCuts; }
  bool isKilledByMissingLayerPair() const { return status_ == Status::killedByMissingLayerPair; }
  bool isKilled() const { return isKilledByCuts() || isKilledByMissingLayerPair(); }
  bool isValidStart() const { return validStart_; }

  // methods to get the vector of inner triplets
  std::vector<Triplet>& innerTriplets();
  std::vector<Triplet> const& innerTripletsView() const;
  size_t innerNeighborIndex(size_t i) const;
  // method to get the number of Triplets
  size_t numInnerTriplets() const;
  // method to get the inner layer ID of the Triplets
  layer_type innerTripletsInnerLayerId() const { return innerTripletsInnerLayerId_; }

private:
  const SimPixelTrack& simPixelTrack_;        // reference to the parent SimPixelTrack
  size_t innerHit_, outerHit_;                // indices of the two hits
  Status status_;                             // status of the doublet
  bool validStart_{false};                    // doublet passes cuts for starting Ntuplets
  int_type numSkippedLayers_;                 // number of layers skipped by the Doublet
  layer_type layerPairId_;                    // ID of the layer pair as defined in the reconstruction for the doublets
  std::vector<Triplet> innerTriplets_;        // indices of inner triplets and its status
  layer_type innerTripletsInnerLayerId_{99};  // layer ID of the inner RecHit of the triplets
};

#endif
