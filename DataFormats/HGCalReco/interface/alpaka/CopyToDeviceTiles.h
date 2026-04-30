#ifndef DataFormats_HGCalReco_interface_alpaka_CopyToDeviceTiles_h
#define DataFormats_HGCalReco_interface_alpaka_CopyToDeviceTiles_h

#include "DataFormats/HGCalReco/interface/TilesHost.h"
#include "DataFormats/HGCalReco/interface/alpaka/TilesDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include <alpaka/alpaka.hpp>

namespace cms::alpakatools {

  // Specialization for copying LayerTiles from host to device
  template <typename T>
  struct CopyToDevice<ticl::LayerTilesHost<T>> {
    template <typename TQueue>
    static auto copyAsync(TQueue& queue, ticl::LayerTilesHost<T> const& src) {
      using TDevice = alpaka::Dev<TQueue>;

      // Create device-side LayerTiles
      // The device type uses ALPAKA_ACCELERATOR_NAMESPACE which expands to the correct device
      using DeviceLayerTiles = ticl::LayerTiles<T, TDevice>;

      // For now, we need to reconstruct the tiles on device
      // This is a placeholder - actual implementation depends on internal structure
      // TODO: Implement proper deep copy of association map data

      DeviceLayerTiles dst(queue, /* size */ 0); // TODO: get actual size

      // Copy internal data structures
      // This will need to copy the association map data
      // alpaka::memcpy(queue, dst.data(), src.data(), size);

      return dst;
    }
  };

  // Specialization for copying full Tiles structure (array of LayerTiles)
  template <typename LayerTilesHost, std::size_t N>
  struct CopyToDevice<ticl::Tiles<LayerTilesHost, N>> {
    template <typename TQueue>
    static auto copyAsync(TQueue& queue, ticl::Tiles<LayerTilesHost, N> const& src) {
      using TDevice = alpaka::Dev<TQueue>;
      using T = typename LayerTilesHost::TilesType;

      // Device-side tiles type
      using DeviceLayerTiles = ticl::LayerTiles<T, TDevice>;
      using DeviceTiles = ticl::Tiles<DeviceLayerTiles, N>;

      // For each layer, copy the LayerTiles to device
      // This creates device tiles with uninitialized data
      DeviceTiles dst(edm::Uninitialized{});

      // Copy each layer's tiles
      // TODO: Implement proper copy for each layer
      // for (size_t i = 0; i < N; ++i) {
      //   // Copy src[i] to dst[i]
      //   // This requires copying the association map data
      // }

      return dst;
    }
  };

}  // namespace cms::alpakatools

#endif
