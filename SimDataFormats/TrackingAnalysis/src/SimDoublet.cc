#include "SimDataFormats/TrackingAnalysis/interface/SimDoublet.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimTriplet.h"

namespace simpixeltracks {

  // Function that determines the number of skipped layers for a given pair of RecHits.
  SimPixelTrack::int_type getNumSkippedLayers(
      std::pair<SimPixelTrack::layer_type, SimPixelTrack::layer_type> const& layerIds,
      unsigned int const inner_detId,
      unsigned int const outer_detId,
      const TrackerTopology* trackerTopology) {
    // Possibility 0: invalid case (outer layer is not the outer one), set to -1 immediately
    if (layerIds.first >= layerIds.second) {
      return -1;
    }

    // get the detector Ids of the two RecHits
    DetId innerDetId(inner_detId);
    DetId outerDetId(outer_detId);

    // determine where the RecHits are
    bool innerInBarrel = (innerDetId.subdetId() == PixelSubdetector::PixelBarrel);
    bool outerInBarrel = (outerDetId.subdetId() == PixelSubdetector::PixelBarrel);
    bool innerInOTExtension = (innerDetId.subdetId() == SiStripSubdetector::TOB);
    bool outerInOTExtension = (outerDetId.subdetId() == SiStripSubdetector::TOB);
    bool innerInBackward = (!innerInBarrel) && (!innerInOTExtension) && (trackerTopology->pxfSide(innerDetId) == 1);
    bool outerInBackward = (!outerInBarrel) && (!outerInOTExtension) && (trackerTopology->pxfSide(outerDetId) == 1);
    bool innerInForward = (!innerInBarrel) && (!innerInOTExtension) && (trackerTopology->pxfSide(innerDetId) == 2);
    bool outerInForward = (!outerInBarrel) && (!outerInOTExtension) && (trackerTopology->pxfSide(outerDetId) == 2);

    // Possibility 1: both RecHits lie in the same detector part (barrel, forward, backward, extension)
    if ((innerInBarrel && outerInBarrel) || (innerInForward && outerInForward) ||
        (innerInBackward && outerInBackward) || (innerInOTExtension && outerInOTExtension)) {
      return (layerIds.second - layerIds.first - 1);
    }
    // Possibility 2: the inner RecHit is in the barrel while the outer is in extension
    else if (innerInBarrel && outerInOTExtension) {
      return (trackerTopology->tobLayer(outerDetId) - trackerTopology->pxbLayer(innerDetId) + 3);
    }
    // Possibility 3: the inner RecHit is in the encaps while the outer is in extension
    else if (outerInOTExtension) {
      return (trackerTopology->tobLayer(outerDetId) - 1);
    }
    // Possibility 4: the inner RecHit is in the barrel while the outer is in either forward or backward
    else if (innerInBarrel) {
      return (trackerTopology->pxfDisk(outerDetId) - 1);
    }
    // Possibility 5: invalid case (one is forward and the other in backward), set to -1
    else {
      return -1;
    }
  }

  // Function that, for a pair of two layers, gives a unique pair Id (innerLayerId * 100 + outerLayerId).
  SimPixelTrack::layer_type getLayerPairId(
      std::pair<SimPixelTrack::layer_type, SimPixelTrack::layer_type> const& layerIds) {
    // calculate the unique layer pair Id as (innerLayerId * 100 + outerLayerId)
    return (layerIds.first * 100 + layerIds.second);
  }
}  // end namespace simpixeltracks

// ------------------------------------------------------------------------------------------------------
// SimPixelTrack::Doublet class member functions
// ------------------------------------------------------------------------------------------------------

// constructor
SimPixelTrack::Doublet::Doublet(SimPixelTrack const& simPixelTrack,
                                size_t const innerIndex,
                                size_t const outerIndex,
                                const TrackerTopology* trackerTopology,
                                std::vector<size_t> const& innerNeighborsIndices)
    : simPixelTrack_(simPixelTrack),
      innerHit_(innerIndex),
      outerHit_(outerIndex),
      // moduleIds_(std::make_pair(simPixelTrack.moduleIds(innerIndex), simPixelTrack.moduleIds(outerIndex))),
      // globalPositions_(
      //     std::make_pair(simPixelTrack.globalPositions(innerIndex), simPixelTrack.globalPositions(outerIndex))),
      // layerIds_(std::make_pair(simPixelTrack.layerIds(innerIndex), simPixelTrack.layerIds(outerIndex))),
      // clusterYSizes_(std::make_pair(simPixelTrack.clusterYSizes(innerIndex), simPixelTrack.clusterYSizes(outerIndex))),
      status_(SimPixelTrack::Doublet::Status::undef) {
  auto layerIds = std::make_pair(simPixelTrack.layerIds(innerIndex), simPixelTrack.layerIds(outerIndex));

  // determine number of skipped layers
  numSkippedLayers_ = simpixeltracks::getNumSkippedLayers(
      layerIds, simPixelTrack.detIds(innerIndex), simPixelTrack.detIds(outerIndex), trackerTopology);

  // determine Id of the layer pair
  layerPairId_ = simpixeltracks::getLayerPairId(layerIds);

  // fill the inner Triplets
  for (size_t const index : innerNeighborsIndices) {
    innerTriplets_.emplace_back(
        SimPixelTrack::Triplet(index, simPixelTrack.getSimDoublet(index).numInnerTriplets()));
  }

  // if there are Triplets, get their inner layerId
  if (!innerNeighborsIndices.empty()) {
    size_t index = innerNeighborsIndices.at(0);
    innerTripletsInnerLayerId_ = simPixelTrack.getSimDoublet(index).innerLayerId();
  }
}

// destructor
SimPixelTrack::Doublet::~Doublet() = default;

std::vector<SimPixelTrack::Triplet>& SimPixelTrack::Doublet::innerTriplets() { return innerTriplets_; }
std::vector<SimPixelTrack::Triplet> const& SimPixelTrack::Doublet::innerTripletsView() const { return innerTriplets_; }
size_t SimPixelTrack::Doublet::innerNeighborIndex(size_t i) const { return innerTriplets_.at(i).innerDoubletIndex(); }
size_t SimPixelTrack::Doublet::numInnerTriplets() const { return innerTriplets_.size(); }
