#pragma once

#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TilesHost.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PatternRecognitionAlgoBase {
  protected:
  public:
    PatternRecognitionAlgoBase(const edm::ParameterSet& conf) {}
    virtual ~PatternRecognitionAlgoBase() = default;

    // Tiles view type - array of LayerTileViews
    using TilesView = ::ticl::TICLLayerTilesHost::TilesType::View;
    using TilesConstViewArray = std::array<TilesView, ::ticl::TileConstants::nLayers>;

    virtual void makeTracksters(Queue& queue,
                                const HGCalSoAClustersDeviceCollection& layerClusters,
                                const TilesConstViewArray& tiles,
                                std::vector<ticl::Trackster>& result) = 0;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
