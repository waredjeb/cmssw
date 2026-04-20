
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/FillAssociator.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include <alpaka/alpaka.hpp>
#include <array>
#include <xtd/xtd.h>
#include <concepts>
#include <cstdint>
#include <span>

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  template <concepts::LayerTile T>
  using LayerTilesDevice = ticl::LayerTiles<T, Device>;

  using TICLLayerTilesHost = Tiles<LayerTilesDevice<TileConstants>, TileConstants::nLayers>;
  using TICLTracksterTilesHost = Tiles<LayerTilesDevice<TileConstants>, TileConstants::iterations>;
  using TICLLayerTilesHFNoseHost = Tiles<LayerTilesDevice<TileConstantsHFNose>, TileConstantsHFNose::nLayers>;
  using TICLTracksterTilesHFNoseHost = Tiles<LayerTilesDevice<TileConstantsHFNose>, TileConstantsHFNose::iterations>;
  using TICLLayerTilesBarrelHost = Tiles<LayerTilesDevice<TileConstantsBarrel>, TileConstantsBarrel::nLayers>;
  using TICLTracksterTilesBarrelHost = Tiles<LayerTilesDevice<TileConstantsBarrel>, TileConstantsBarrel::iterations>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
