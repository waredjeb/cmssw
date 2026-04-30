
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/HGCalReco/interface/Tiles.h"
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

  using TICLLayerTilesDevice = ::ticl::Tiles<::ticl::TileConstants, ::ticl::TileConstants::nLayers, Device>;
  using TICLTracksterTilesDevice = ::ticl::Tiles<::ticl::TileConstants, ::ticl::TileConstants::iterations, Device>;
  using TICLLayerTilesHFNoseDevice = ::ticl::Tiles<::ticl::TileConstantsHFNose, ::ticl::TileConstantsHFNose::nLayers, Device>;
  using TICLTracksterTilesHFNoseDevice = ::ticl::Tiles<::ticl::TileConstantsHFNose, ::ticl::TileConstantsHFNose::iterations, Device>;
  using TICLLayerTilesBarrelDevice = ::ticl::Tiles<::ticl::TileConstantsBarrel, ::ticl::TileConstantsBarrel::nLayers, Device>;
  using TICLTracksterTilesBarrelDevice = ::ticl::Tiles<::ticl::TileConstantsBarrel, ::ticl::TileConstantsBarrel::iterations, Device>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
