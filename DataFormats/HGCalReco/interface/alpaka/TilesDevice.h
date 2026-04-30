
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

  using TICLLayerTilesHost = Tiles<TileConstants, TileConstants::nLayers, Device>;
  using TICLTracksterTilesHost = Tiles<TileConstants, TileConstants::iterations, Device>;
  using TICLLayerTilesHFNoseHost = Tiles<TileConstantsHFNose, TileConstantsHFNose::nLayers, Device>;
  using TICLTracksterTilesHFNoseHost = Tiles<TileConstantsHFNose, TileConstantsHFNose::iterations, Device>;
  using TICLLayerTilesBarrelHost = Tiles<TileConstantsBarrel, TileConstantsBarrel::nLayers, Device>;
  using TICLTracksterTilesBarrelHost = Tiles<TileConstantsBarrel, TileConstantsBarrel::iterations, Device>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
