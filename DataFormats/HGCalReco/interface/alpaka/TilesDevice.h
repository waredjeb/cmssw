
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/HGCalReco/interface/Tiles.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  // Device-side layer tiles (using the backend-specific Device type)
  template <concepts::LayerTile T>
  using LayerTilesDevice = ::ticl::LayerTiles<T, Device>;

  // Device-side tiles collections
  using TICLLayerTilesDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstants>, ::ticl::TileConstants::nLayers>;
  using TICLTracksterTilesDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstants>, ::ticl::TileConstants::iterations>;
  using TICLLayerTilesHFNoseDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstantsHFNose>, ::ticl::TileConstantsHFNose::nLayers>;
  using TICLTracksterTilesHFNoseDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstantsHFNose>, ::ticl::TileConstantsHFNose::iterations>;
  using TICLLayerTilesBarrelDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstantsBarrel>, ::ticl::TileConstantsBarrel::iterations>;
  using TICLTracksterTilesBarrelDevice = ::ticl::Tiles<LayerTilesDevice<::ticl::TileConstantsBarrel>, ::ticl::TileConstantsBarrel::iterations>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
