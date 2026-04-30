
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/HGCalReco/interface/Tiles.h"
#include <alpaka/alpaka.hpp>

namespace ticl {

  using TICLLayerTilesHost = Tiles<TileConstants, TileConstants::nLayers, alpaka::DevCpu>;
  using TICLTracksterTilesHost = Tiles<TileConstants, TileConstants::iterations, alpaka::DevCpu>;
  using TICLLayerTilesHFNoseHost = Tiles<TileConstantsHFNose, TileConstantsHFNose::nLayers, alpaka::DevCpu>;
  using TICLTracksterTilesHFNoseHost = Tiles<TileConstantsHFNose, TileConstantsHFNose::iterations, alpaka::DevCpu>;
  using TICLLayerTilesBarrelHost = Tiles<TileConstantsBarrel, TileConstantsBarrel::nLayers, alpaka::DevCpu>;
  using TICLTracksterTilesBarrelHost = Tiles<TileConstantsBarrel, TileConstantsBarrel::iterations, alpaka::DevCpu>;
  using TICLTracksterLinkingTilesHost = Tiles<TileConstants, 2, alpaka::DevCpu>;

}  // namespace ticl
