#include "SimDataFormats/TrackingAnalysis/interface/SimPixelTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimDoublet.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimTriplet.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimNtuplet.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

// default contructor
SimPixelTrack::SimPixelTrack() = default;

// constructor
SimPixelTrack::SimPixelTrack(TrackingParticleRef const trackingParticleRef, reco::BeamSpot const& beamSpot)
    : trackingParticleRef_(trackingParticleRef), beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}
SimPixelTrack::SimPixelTrack(reco::TrackBaseRef const trackRef, reco::BeamSpot const& beamSpot)
    : trackRef_(trackRef), beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}
SimPixelTrack::SimPixelTrack(reco::BeamSpot const& beamSpot)
    : beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}

// destructor
SimPixelTrack::~SimPixelTrack() = default;

// method to add a RecHit to the SimPixelTrack
void SimPixelTrack::addRecHit(TrackingRecHit const& recHit,
                              layer_type const layerId,
                              int_type const clusterYSize,
                              unsigned int const detId,
                              int const moduleId) {
  recHitsAreSorted_ = false;  // set sorted-bool to false again

  // check if the layerId is not present in the layerIdVector yet
  if (std::find(hits_.layerIds.begin(), hits_.layerIds.end(), layerId) == hits_.layerIds.end()) {
    // if it does not exist, increment number of layers
    numLayers_++;
  }

  // add detId, the corrected hit position, layerId and clusterSize to respective vectors
  hits_.detIds.push_back(detId);
  hits_.moduleIds.push_back(moduleId);
  hits_.globalPositions.push_back(recHit.globalPosition() - beamSpotPosition_);
  hits_.layerIds.push_back(layerId);
  hits_.clusterYSizes.push_back(clusterYSize);
}

// method to sort the RecHits according to the position relative to the TP vertex
void SimPixelTrack::sortRecHits() {
  auto vertex = trackingParticleRef_->vertex();
  sortRecHits(vertex.x(), vertex.y(), vertex.z());
}
// method to sort the RecHits according to the position relative to a given reference
void SimPixelTrack::sortRecHits(float const x, float const y, float const z) {
  // get the production vertex of the TrackingParticle (corrected for beamspot)
  const GlobalVector vertex(x - beamSpotPosition_.x(), y - beamSpotPosition_.y(), z - beamSpotPosition_.z());

  // get the vector of squared magnitudes of the global RecHit positions relative to vertex
  std::vector<double> recHitMag2;
  recHitMag2.reserve(hits_.layerIds.size());
  for (const auto& globalPosition : hits_.globalPositions) {
    // relative RecHit position with respect to the production vertex
    Global3DPoint relativePosition = globalPosition - vertex;
    recHitMag2.push_back(relativePosition.mag2());
  }

  // find the permutation vector that sorts the magnitudes
  std::vector<std::size_t> sortedPerm(recHitMag2.size());
  std::iota(sortedPerm.begin(), sortedPerm.end(), 0);
  std::sort(sortedPerm.begin(), sortedPerm.end(), [&](std::size_t i, std::size_t j) {
    return (recHitMag2[i] < recHitMag2[j]);
  });

  // create the sorted vectors
  std::vector<unsigned int> sorted_detIdVector;
  std::vector<int> sorted_moduleIdVector;
  std::vector<GlobalPoint> sorted_globalPositionVector;
  std::vector<layer_type> sorted_layerIdVector;
  std::vector<int_type> sorted_clusterYSizeVector;
  sorted_detIdVector.reserve(sortedPerm.size());
  sorted_moduleIdVector.reserve(sortedPerm.size());
  sorted_globalPositionVector.reserve(sortedPerm.size());
  sorted_layerIdVector.reserve(sortedPerm.size());
  sorted_clusterYSizeVector.reserve(sortedPerm.size());
  for (size_t i : sortedPerm) {
    sorted_detIdVector.push_back(hits_.detIds[i]);
    sorted_moduleIdVector.push_back(hits_.moduleIds[i]);
    sorted_globalPositionVector.push_back(hits_.globalPositions[i]);
    sorted_layerIdVector.push_back(hits_.layerIds[i]);
    sorted_clusterYSizeVector.push_back(hits_.clusterYSizes[i]);
  }

  // swap them with the class member
  hits_.detIds.swap(sorted_detIdVector);
  hits_.moduleIds.swap(sorted_moduleIdVector);
  hits_.globalPositions.swap(sorted_globalPositionVector);
  hits_.layerIds.swap(sorted_layerIdVector);
  hits_.clusterYSizes.swap(sorted_clusterYSizeVector);

  // set sorted bool to true
  recHitsAreSorted_ = true;
}

// method to get fishbone alignments
std::vector<SimPixelTrack::Fishbone> SimPixelTrack::fishboneScores() const {
  // confirm that the RecHits are sorted
  assert(recHitsAreSorted_);

  std::vector<Fishbone> fishbones{};

  if (numRecHits() < 3) {
    return fishbones;
  }

  // loop over outer hits
  for (size_t o{0}; o < numRecHits(); o++) {
    auto outerLayerId = layerIds(o);
    auto outerPos = globalPositions(o);
    auto nInnerDoublets = innerDoubletsOfRecHit_.at(o).size();

    // loop over first inner doublets of that outer hit
    for (size_t d1{0}; d1 < nInnerDoublets; d1++) {
      auto innerPos1 = getSimDoublet(d1).innerGlobalPos();
      double x1 = (innerPos1.x() - outerPos.x());
      double y1 = (innerPos1.y() - outerPos.y());
      double z1 = (innerPos1.z() - outerPos.z());
      double n1 = x1 * x1 + y1 * y1 + z1 * z1;
      // loop over first inner doublets of that outer hit
      for (size_t d2{d1 + 1}; d2 < nInnerDoublets; d2++) {
        auto innerPos2 = getSimDoublet(d2).innerGlobalPos();
        double x2 = (innerPos2.x() - outerPos.x());
        double y2 = (innerPos2.y() - outerPos.y());
        double z2 = (innerPos2.z() - outerPos.z());
        double n2 = x2 * x2 + y2 * y2 + z2 * z2;
        auto cos12 = x1 * x2 + y1 * y2 + z1 * z2;

        auto fishboneCut = cos12 * cos12 / (n1 * n2);

        fishbones.emplace_back(Fishbone(outerLayerId, fishboneCut));
      }
    }
  }
  return fishbones;
}

// method to produce the true doublets
void SimPixelTrack::buildSimDoublets(const TrackerTopology* trackerTopology) const {
  // confirm that the RecHits are sorted
  assert(recHitsAreSorted_);

  // check if there are at least two hits
  if (numRecHits() < 2) {
    return;
  }

  // resize vector innerDoubletsOfRecHit_ to actual number of RecHits
  innerDoubletsOfRecHit_.resize(numRecHits());

  // updatable current number of doublets
  size_t nDoublets{0};

  // loop over the RecHits/layer Ids
  for (size_t i = 0; i < hits_.layerIds.size(); i++) {
    layer_type innerLayerId = hits_.layerIds[i];
    layer_type outerLayerId{};
    size_t outerLayerStart{hits_.layerIds.size()};

    // find the next layer Id + at which hit this layer starts
    for (size_t j = i + 1; j < hits_.layerIds.size(); j++) {
      if (innerLayerId != hits_.layerIds[j]) {
        outerLayerId = hits_.layerIds[j];
        outerLayerStart = j;
        break;
      }
    }

    // build the doublets of the inner hit i with all outer hits j in the layer outerLayerId
    for (size_t j = outerLayerStart; j < hits_.layerIds.size(); j++) {
      // break if the hit doesn't belong to the outer layer anymore
      if (outerLayerId != hits_.layerIds[j]) {
        break;
      }

      // create and append new doublet
      doublets_.emplace_back(SimPixelTrack::Doublet(*this, i, j, trackerTopology, innerDoubletsOfRecHit_.at(i)));

      // save the index of the new doublet in the outer RecHit's innerDoubletsOfRecHit_
      innerDoubletsOfRecHit_.at(j).push_back(nDoublets);

      // update the number of doublets
      nDoublets++;
    }
  }  // end loop over the RecHits/layer Ids
}

// function to recursively build the Ntuplets from a given starting doublet
// (the building starts from the outside and ends inside)
// at each addition of a SimDoublet, a new SimNtuplet is stored
void SimPixelTrack::buildSimNtuplets(Doublet const& doublet,
                                     std::vector<bool> const& quadruplets,
                                     size_t numSimDoublets,
                                     layer_type const lastLayerId,
                                     status_type const status,
                                     int_type const numSkippedLayers,
                                     size_t const minNumDoubletsToPass) const {
  // update the number of SimDoublets once before looping over the actual neighbors to be added
  numSimDoublets++;

  // loop over the inner neighboring doublets of the current doublet
  for (size_t i{0}; auto const& triplet : doublet.innerTripletsView()) {
    // get the inner neighboring doublet and the status of this connection
    auto const& neighborDoublet = doublets_.at(triplet.innerDoubletIndex());

    // update the status of the current SimNtuplet by adding the information from the new doublet
    status_type updatedStatus = SimPixelTrack::Ntuplet::updateStatus(
        status,                                           // current status
        neighborDoublet.isUndef(),                        // doublet has undefined cuts
        neighborDoublet.isKilledByMissingLayerPair(),     // doublet is not built due to missing layer pair
        neighborDoublet.isKilledByCuts(),                 // doublet is killed by cuts
        triplet.isUndef(),                                // doublet connection has undefined cuts
        triplet.isKilled(),                               // doublet connection is killed by cuts
        (numSimDoublets > 2) ? quadruplets.at(i) : false  // triplet connection is killed by cuts
    );

    // update number of skipped layers
    int_type updatedNumSkippedLayers = numSkippedLayers + doublet.numSkippedLayers();

    // add the current state as a new SimNtuplet to the collection
    ntuplets_.emplace_back(SimPixelTrack::Ntuplet(numSimDoublets,
                                                  updatedStatus,
                                                  neighborDoublet.innerLayerId(),
                                                  neighborDoublet.outerLayerId(),
                                                  lastLayerId,
                                                  updatedNumSkippedLayers));

    // change the status "TooShort" of the newly created SimNtuplet if it is indeed to short
    if (numSimDoublets < minNumDoubletsToPass) {
      ntuplets_.back().setTooShort();
    }

    // change the status "invalidStart" of the newly created SimNtuplet if this is indeed the case
    if (!(neighborDoublet.isValidStart())) {
      ntuplets_.back().setInvalidStart();
    }

    // check if the new SimNtuplet qualifies as longest SimNtuplet
    // A) if it's the first Ntuplet or longer than the current longest,
    //    it becomes automatically the longest
    // B) otherwise:
    //     - it needs to be at least as long as the current longest
    //     - and it needs to get farther in the reconstruction chain:
    //        1. Ntuplet is long enough
    //        2. no missing layer pairs
    //        3. all doublets survive
    //        4. all doublet connections survive
    //        5. all triplet connections survive
    //        6. first doublet from starting layer pair
    if ((!longestNtupletIndex_) || (numSimDoublets > ntuplets_.at(*longestNtupletIndex_).numDoublets())) {
      // case A)
      longestNtupletIndex_ = ntuplets_.size() - 1;
    } else if ((numSimDoublets == ntuplets_.at(*longestNtupletIndex_).numDoublets()) &&  // is at least as long
               ntuplets_.back().getsFartherInRecoChainThanReference(
                   ntuplets_.at(*longestNtupletIndex_))) {  // get farther in reconstruction
      // case B)
      longestNtupletIndex_ = ntuplets_.size() - 1;
    }

    // check if the new SimNtuplet qualifies as longest SimNtuplet alive
    if (ntuplets_.back().isAlive()) {      // obviously, it has to be alive
      if ((!longestAliveNtupletIndex_) ||  // it's the first SimNtuplet alive or
          ((numSimDoublets >= ntuplets_.at(*longestAliveNtupletIndex_).numDoublets()) &&  // is at least as long and
           (ntuplets_.back().firstLayerId() <=
            ntuplets_.at(*longestAliveNtupletIndex_).firstLayerId()))  // is at least as inside
      ) {
        longestAliveNtupletIndex_ = ntuplets_.size() - 1;
      }
    }

    // check if the new SimNtuplet qualifies as best SimNtuplet yet
    // A) if it's the first Ntuplet or farther in the reco chain than the current best,
    //    it becomes automatically the best
    // B) otherwise:
    //     - it needs to be at least as long as the current longest
    //     - and it needs to get at least as far in the reconstruction chain:
    //        1. Ntuplet is long enough
    //        2. no missing layer pairs
    //        3. all doublets survive
    //        4. all doublet connections survive
    //        5. all triplet connections survive
    //        6. first doublet from starting layer pair
    if ((!bestNtupletIndex_) ||
        ntuplets_.back().getsFartherInRecoChainThanReference(ntuplets_.at(*bestNtupletIndex_))) {
      // case A)
      bestNtupletIndex_ = ntuplets_.size() - 1;
    } else if ((numSimDoublets >= ntuplets_.at(*bestNtupletIndex_).numDoublets()) &&  // is at least as long
               ntuplets_.back().getsAsFarInRecoChainAsReference(
                   ntuplets_.at(*bestNtupletIndex_))) {  // get as far in reconstruction
      // case B)
      bestNtupletIndex_ = ntuplets_.size() - 1;
    }

    // call this function recursively
    // this will get the further neighboring doublets and build the next Ntuplet
    buildSimNtuplets(neighborDoublet,
                     triplet.quadruplets(),
                     numSimDoublets,
                     lastLayerId,
                     updatedStatus,
                     updatedNumSkippedLayers,
                     minNumDoubletsToPass);
    i++;
  }
}

// method to produce the SimNtuplets
// (collection of all possible Ntuplets you can build from the SimDoublets)
void SimPixelTrack::buildSimNtuplets(size_t const minNumDoubletsToPass) const {
  // clear the Ntuplet collection and reset longest Ntuplet indices
  ntuplets_.clear();
  longestNtupletIndex_.reset();
  longestAliveNtupletIndex_.reset();
  bestNtupletIndex_.reset();

  // check if there are at least two doublets
  if (numDoublets() < 2) {
    return;
  }

  // loop over all SimDoublets, using them as starting points for building Ntuplets
  for (auto const& doublet : doublets_) {
    // intialize status according to the doublet properties
    status_type status = SimPixelTrack::Ntuplet::updateStatus(
        0,                                     // current status to be updated
        doublet.isUndef(),                     // doublet has undefined cuts
        doublet.isKilledByMissingLayerPair(),  // doublet is not built due to missing layer pair
        doublet.isKilledByCuts(),              // doublet is killed by cuts
        false,                                 // doublet connection has undefined cuts
        false,                                 // doublet connection is killed by cuts
        false                                  // triplet connection is killed by cuts
    );
    // initialize number of skipped layers
    int_type numSkippedLayers = doublet.numSkippedLayers();
    // build the Ntuplets recursively
    buildSimNtuplets(doublet, {}, 1, doublet.outerLayerId(), status, numSkippedLayers, minNumDoubletsToPass);
  }
}

size_t SimPixelTrack::numDoublets() const { return doublets_.size(); }

SimPixelTrack::Doublet const& SimPixelTrack::getSimDoublet(size_t const index) const { return doublets_.at(index); }

std::vector<SimPixelTrack::Doublet>& SimPixelTrack::getSimDoublets() const { return doublets_; }
std::vector<SimPixelTrack::Doublet>& SimPixelTrack::buildAndGetSimDoublets(
    const TrackerTopology* trackerTopology) const {
  buildSimDoublets(trackerTopology);
  return doublets_;
}

std::vector<SimPixelTrack::Ntuplet>& SimPixelTrack::getSimNtuplets() const { return ntuplets_; };
std::vector<SimPixelTrack::Ntuplet>& SimPixelTrack::buildAndGetSimNtuplets(size_t const minNumDoubletsToPass = 0) const {
  buildSimNtuplets(minNumDoubletsToPass);
  return ntuplets_;
};

bool SimPixelTrack::hasSimNtuplet() const { return longestNtupletIndex_.has_value(); }
bool SimPixelTrack::hasAliveSimNtuplet() const { return longestAliveNtupletIndex_.has_value(); }

SimPixelTrack::Ntuplet const& SimPixelTrack::longestSimNtuplet() const { return ntuplets_.at(*longestNtupletIndex_); }
SimPixelTrack::Ntuplet const& SimPixelTrack::longestAliveSimNtuplet() const {
  return ntuplets_.at(*longestAliveNtupletIndex_);
}
SimPixelTrack::Ntuplet const& SimPixelTrack::bestSimNtuplet() const { return ntuplets_.at(*bestNtupletIndex_); }

void SimPixelTrack::clearMutables() const {
  doublets_.clear();
  ntuplets_.clear();
  longestNtupletIndex_.reset();
  longestAliveNtupletIndex_.reset();
  bestNtupletIndex_.reset();
}
