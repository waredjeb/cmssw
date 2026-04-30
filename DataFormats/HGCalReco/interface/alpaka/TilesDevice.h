
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

  using TICLLayerTilesDevice = Tiles<TileConstants, TileConstants::nLayers, Device>;
  using TICLTracksterTilesDevice = Tiles<TileConstants, TileConstants::iterations, Device>;
  using TICLLayerTilesHFNoseDevice = Tiles<TileConstantsHFNose, TileConstantsHFNose::nLayers, Device>;
  using TICLTracksterTilesHFNoseDevice = Tiles<TileConstantsHFNose, TileConstantsHFNose::iterations, Device>;
  using TICLLayerTilesBarrelDevice = Tiles<TileConstantsBarrel, TileConstantsBarrel::nLayers, Device>;
  using TICLTracksterTilesBarrelDevice = Tiles<TileConstantsBarrel, TileConstantsBarrel::iterations, Device>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
