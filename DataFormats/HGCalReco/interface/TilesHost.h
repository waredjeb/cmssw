
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/HGCalReco/interface/Tiles.h"
#include <alpaka/alpaka.hpp>

namespace ticl {

  template <typename T>
  using LayerTilesHost = ticl::LayerTiles<T, alpaka::DevCpu>;

  using TICLLayerTilesHost = Tiles<LayerTilesHost<TileConstants>, TileConstants::nLayers>;
  using TICLTracksterTilesHost = Tiles<LayerTilesHost<TileConstants>, TileConstants::iterations>;
  using TICLLayerTilesHFNoseHost = Tiles<LayerTilesHost<TileConstantsHFNose>, TileConstantsHFNose::nLayers>;
  using TICLTracksterTilesHFNoseHost = Tiles<LayerTilesHost<TileConstantsHFNose>, TileConstantsHFNose::iterations>;
  using TICLLayerTilesBarrelHost = Tiles<LayerTilesHost<TileConstantsBarrel>, TileConstantsBarrel::nLayers>;
  using TICLTracksterTilesBarrelHost = Tiles<LayerTilesHost<TileConstantsBarrel>, TileConstantsBarrel::iterations>;
  using TICLTracksterLinkingTilesHost = Tiles<LayerTilesHost<TileConstants>, 2>;

}  // namespace ticl
