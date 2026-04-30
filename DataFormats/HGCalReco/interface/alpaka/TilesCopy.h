#ifndef DataFormats_HGCalReco_alpaka_TilesCopy_h
#define DataFormats_HGCalReco_alpaka_TilesCopy_h

#include "DataFormats/HGCalReco/interface/TilesHost.h"
#include "DataFormats/HGCalReco/interface/alpaka/TilesDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  // Copy tiles from host to device
  // This function creates device-side tiles and copies the data asynchronously
  template <typename TQueue>
  auto copyTilesToDevice(TQueue& queue, const ::ticl::TICLLayerTilesHost& hostTiles) {
    // Create device tiles with same structure
    TICLLayerTilesDevice deviceTiles(edm::Uninitialized{});

    // Copy each layer's tiles from host to device
    // Note: The actual copy implementation depends on the internal structure
    // of PortableCollection. For now, this is a placeholder.
    //
    // TODO: Implement proper copy of PortableCollection data
    // This will likely involve:
    // 1. For each layer in hostTiles
    // 2. Get the PortableCollection from hostTiles[i]
    // 3. Create corresponding device PortableCollection
    // 4. Copy data using alpaka::memcpy
    //
    // The challenge is that PortableCollection doesn't expose a direct copy constructor
    // that takes host collection and creates device collection.

    return deviceTiles;
  }

  // View type that can be passed to kernels
  template <typename LayerTilesDevice, std::size_t N>
  using TilesView = std::array<typename LayerTilesDevice::View, N>;

  // Get a view of device tiles that can be used in kernels
  template <typename LayerTilesDevice, std::size_t N>
  auto getTilesView(const ::ticl::Tiles<LayerTilesDevice, N>& tiles) {
    return tiles.view();
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl

#endif
