#ifndef SimDataFormats_TrackingAnalysis_SimPixelTrack_h
#define SimDataFormats_TrackingAnalysis_SimPixelTrack_h

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitFwd.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/** @brief Semi-Monte Carlo truth information used for pixel-tracking opimization.
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
 * @date January 2025
 */
class SimPixelTrack {
public:
  // types for SimPixelTrack properties (also used for CA cuts in the SimPixelTracksAnalyzer)
  using float_type = double;
  using int_type = int16_t;
  using layer_type = uint16_t;
  using status_type = uint8_t;

  inline static layer_type kInvalidLayerId = std::numeric_limits<layer_type>::max();

  struct Hits {
    // vectors of usual hit properties
    std::vector<unsigned int> detIds;          // detector Ids of the RecHits
    std::vector<int> moduleIds;                // module Ids of the RecHits
    std::vector<GlobalPoint> globalPositions;  // global positions of the RecHits (corrected by beamspot)
    std::vector<layer_type> layerIds;          // layer IDs corresponding to the RecHits
    std::vector<int_type> clusterYSizes;       // cluster sizes (local y) corresponding to the RecHits

    // vectors of vectorHit properties
    std::vector<float_type> dPhiDrs;     // dPhi / dR of the stub, -1 for non-stub hit
    std::vector<float_type> dPhiDrErrs;  // error of dPhi / dR of the stub, -1 for non-stub hit
  };

  struct Fishbone {
    layer_type layerId;  // outer layer Id of the doublets
    float_type score;    // fishbone score
  };

  // SimDoublets, SimTriplets and SimNtuplets
  class Doublet;
  class Triplet;
  class Ntuplet;

  // contructors
  SimPixelTrack();
  SimPixelTrack(TrackingParticleRef const trackingParticleRef, reco::BeamSpot const& beamSpot);
  SimPixelTrack(reco::TrackBaseRef const trackRef, reco::BeamSpot const& beamSpot);
  SimPixelTrack(reco::BeamSpot const& beamSpot);

  // destructor
  ~SimPixelTrack();

  // method to add a RecHit to the SimPixelTrack
  void addRecHit(TrackingRecHit const& recHit,
                 layer_type const layerId,
                 int_type const clusterYSize,
                 unsigned int const detId,
                 int const moduleId);

  // method to get the reference to the TrackingParticle
  TrackingParticleRef trackingParticle() const { return trackingParticleRef_; }
  // method to get the reference to the track
  reco::TrackBaseRef track() const { return trackRef_; }

  // method to get the detector id vector
  std::vector<unsigned int> detIds() const { return hits_.detIds; }
  // method to get the detector id at index i
  unsigned int detIds(size_t const i) const { return hits_.detIds[i]; }

  // method to get the module id vector
  std::vector<int> moduleIds() const { return hits_.moduleIds; }
  // method to get the module id at index i
  int moduleIds(size_t const i) const { return hits_.moduleIds[i]; }

  // method to get the global position vector of the RecHits
  std::vector<GlobalPoint> globalPositions() const { return hits_.globalPositions; }
  // method to get the global position of the RecHit at index i
  GlobalPoint globalPositions(size_t const i) const { return hits_.globalPositions[i]; }

  // method to get the layer id vector
  std::vector<layer_type> layerIds() const { return hits_.layerIds; }
  // method to get the layer id at index i
  layer_type layerIds(size_t const i) const { return hits_.layerIds[i]; }

  // method to get the cluster size vector
  std::vector<int_type> clusterYSizes() const { return hits_.clusterYSizes; }
  // method to get the cluster size at index i
  int_type clusterYSizes(size_t const i) const { return hits_.clusterYSizes[i]; }

  // method to get the beam spot position
  GlobalVector beamSpotPosition() const { return beamSpotPosition_; }

  // method to get the number of layers
  size_t numLayers() const { return numLayers_; }
  // method to get number of RecHits in the SimPixelTrack
  size_t numRecHits() const { return hits_.layerIds.size(); }
  // method to get the number of SimDoublets
  size_t numDoublets() const;

  // method to sort the RecHits according to the position (either a given reference point or the TP vertex)
  void sortRecHits();
  void sortRecHits(float const, float const, float const);

  // method to produce the SimDoublets from the RecHits
  void buildSimDoublets(const TrackerTopology* trackerTopology) const;
  // method to access the SimDoublets
  std::vector<Doublet>& getSimDoublets() const;
  // method to build and access the SimDoublets
  std::vector<Doublet>& buildAndGetSimDoublets(const TrackerTopology* trackerTopology) const;
  // method to access a single SimDoublet
  Doublet const& getSimDoublet(size_t const index) const;

  // method to build the SimNtuplets
  // minNumDoubletsToPass = the number of doublets required for the Ntuplet to not be considered too short
  void buildSimNtuplets(size_t const minNumDoubletsToPass) const;
  // method to access the SimNtuplets
  std::vector<Ntuplet>& getSimNtuplets() const;
  // method to build and access the SimNtuplets in one go
  // minNumDoubletsToPass = the number of doublets required for the Ntuplet to not be considered too short
  std::vector<Ntuplet>& buildAndGetSimNtuplets(size_t const minNumDoubletsToPass) const;

  // method to check if there are SimNtuplets
  bool hasSimNtuplet() const;
  // method to check if there are alive SimNtuplet
  bool hasAliveSimNtuplet() const;

  // method to access the longest SimNtuplet
  Ntuplet const& longestSimNtuplet() const;
  // method to access the longest alive SimNtuplet
  Ntuplet const& longestAliveSimNtuplet() const;
  // method to access the best SimNtuplet
  Ntuplet const& bestSimNtuplet() const;

  // method to get fishbone alignments
  std::vector<Fishbone> fishboneScores() const;

  // method to clear the mutable vectors once you finished using them
  void clearMutables() const;

private:
  // function for recursive building of Ntuplets
  void buildSimNtuplets(Doublet const& doublet,
                        std::vector<bool> const& quadruplets,
                        size_t numSimDoublets,
                        layer_type const lastLayerId,
                        status_type const status,
                        int_type const numSkippedLayers,
                        size_t const minNumDoubletsToPass) const;

  // class members
  TrackingParticleRef trackingParticleRef_;  // reference to the TrackingParticle (if SimPixelTrack is based on a TP)
  reco::TrackBaseRef trackRef_;              // referency to the track (if SimPixelTrack is based on a track)
  Hits hits_;                                // RecHits associated to the TP
  GlobalVector beamSpotPosition_;  // global position of the beam spot (needed to correct the global RecHit position)
  bool recHitsAreSorted_{false};   // true if RecHits were sorted
  size_t numLayers_{0};            // number of layers hit by the TrackingParticle

  // non-persistent, mutable members:
  // vector of true doublets
  mutable std::vector<Doublet> doublets_;
  // vector of true Ntuplets
  mutable std::vector<Ntuplet> ntuplets_;
  // index of the longest SimNtuplet
  mutable std::optional<size_t> longestNtupletIndex_{-1};
  // index of the longest SimNtuplet that survives
  mutable std::optional<size_t> longestAliveNtupletIndex_{-1};
  // index of the SimNtuplet that gets the farthest in the reco chain
  mutable std::optional<size_t> bestNtupletIndex_{-1};
  // vector of length NrecHits that holds for each RecHit references to all the
  // doublets that have this hit as an outer hit
  mutable std::vector<std::vector<size_t>> innerDoubletsOfRecHit_{};
};

// collection of SimPixelTrack
typedef std::vector<SimPixelTrack> SimPixelTrackCollection;

#endif
