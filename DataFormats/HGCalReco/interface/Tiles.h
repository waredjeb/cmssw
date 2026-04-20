
#pragma once

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/Portable/interface/PortableCollection.h"
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

namespace ticl {

  template <concepts::LayerTile T>
  class LayerTilesView {
  private:
    ::ticl::AssociationMapConstView<> associations;

  public:
    constexpr LayerTilesView() = default;
    constexpr LayerTilesView(::ticl::AssociationMapConstView<> assoc_view) : associations{assoc_view} {}

    ALPAKA_FN_ACC auto operator[](std::size_t idx) const { return associations[idx]; }
    ALPAKA_FN_ACC auto operator[](std::size_t idx) { return associations[idx]; }

    ALPAKA_FN_ACC auto contains(std::size_t idx) const { return associations.contains(idx); }
    ALPAKA_FN_ACC auto count(std::size_t idx) const { return associations.count(idx); }

    ALPAKA_FN_ACC auto phiBin(float phi) const {
      const auto normPhi = normalizedPhi(phi);
      const auto r = T::nPhiBins * M_1_PI * 0.5f;
      const auto phiBin = static_cast<int>((normPhi + M_PI) * r);

      return phiBin;
    }
    ALPAKA_FN_ACC auto etaBin(float eta) const {
      constexpr auto etaRange = T::maxEta - T::minEta;
      static_assert(etaRange >= 0.f);
      const auto r = T::nEtaBins / etaRange;
      int etaBin;
      // if constexpr (std::is_same_v<T, ticl::TileConstantsBarrel>)
      //   etaBin = (eta - T::minEta) * r;
      // else
      etaBin = (xtd::abs(eta) - T::minEta) * r;
      etaBin = xtd::clamp(etaBin, 0, T::nEtaBins - 1);
      return etaBin;
    }

    ALPAKA_FN_ACC auto globalBin(int etaBin, int phiBin) const { return phiBin + etaBin * T::nPhiBins; }
    ALPAKA_FN_ACC auto globalBin(float eta, float phi) const { return phiBin(phi) + etaBin(eta) * T::nPhiBins; }

    ALPAKA_FN_ACC auto searchBox(float etaMin, float etaMax, float phiMin, float phiMax) const {
      // if (!std::is_same_v<T, ticl::TileConstantsBarrel>) {
      //   if (etaMin * etaMax < 0) {
      //     return std::array<int, 4>({{0, 0, 0, 0}});
      //   }
      // }
      if (etaMax - etaMin < 0) {
        return std::array<int, 4>({{0, 0, 0, 0}});
      }
      auto etaBinMin = etaBin(etaMin);
      auto etaBinMax = etaBin(etaMax);
      const auto phiBinMin = phiBin(phiMin);
      auto phiBinMax = phiBin(phiMax);
      if (etaMin < 0) {
        std::swap(etaBinMin, etaBinMax);
      }

      if (phiBinMax < phiBinMin) {
        phiBinMax += T::nPhiBins;
      }
      return std::array<int, 4>({{etaBinMin, etaBinMax, phiBinMin, phiBinMax}});
    }
  };

  struct KernelTilesAssociations {
    template <alpaka::concepts::Acc TAcc, concepts::LayerTile T>
    ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                  LayerTilesView<T> tiles,
                                  std::span<const float> etas,
                                  std::span<const float> phis,
                                  uint32_t* associations) const {
      for (auto idx : alpaka::uniformElements(acc, etas.size())) {
        associations[idx] = tiles.globalBin(etas[idx], phis[idx]);
      }
    }
  };

  // TODO: check consistensy of integer signedness
  template <concepts::LayerTile T, typename TDev>
  class LayerTiles {
  public:
    using View = LayerTilesView<T>;
    using TilesType = T;

    explicit LayerTiles() noexcept : m_associations(edm::Uninitialized{}) {}
    explicit LayerTiles(edm::Uninitialized init) noexcept : m_associations(init) {}
    template <typename TQueue>
    LayerTiles(TQueue& queue, std::integral auto nvalues) : m_associations(queue, nvalues, T::nBins) {}
    LayerTiles(std::integral auto nvalues)
      requires std::same_as<TDev, alpaka::DevCpu>
        : m_associations(cms::alpakatools::host(), nvalues, T::nBins) {}

    template <typename TAcc, typename TQueue>
    ALPAKA_FN_HOST void fill(TQueue& queue,
                             std::span<const float> etas,
                             std::span<const float> phis,
                             std::span<const uint32_t> layer_cluster_ids) {
      auto associations = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, etas.size());

      const auto blocksize = 1024u;
      const auto gridsize = cms::alpakatools::divide_up_by(etas.size(), blocksize);
      auto work_division = cms::alpakatools::make_workdiv<TAcc>(gridsize, blocksize);
      alpaka::exec<TAcc>(
          queue, work_division, KernelTilesAssociations{}, this->view(), etas, phis, associations.data());

      associator::fill<TAcc>(
          queue, m_associations.view(), std::span<const uint32_t>(associations.data(), etas.size()), layer_cluster_ids);
    }

    ALPAKA_FN_HOST auto view() const { return View(m_associations.view()); }
    ALPAKA_FN_HOST auto view() { return View(m_associations.view()); }

    // TODO: do we want to be able to clear them?
    // void clear(Queue& queue);

  private:
    PortableCollection<TDev, ::ticl::AssociationMap<>> m_associations;
  };

  template <typename LayerTiles, std::size_t N>
  class Tiles {
  public:
    using View = LayerTiles::View;
    using LayerTilesType = LayerTiles;
    using TilesType = LayerTiles::TilesType;

    explicit Tiles(edm::Uninitialized init) noexcept {
      for (auto dim = 0u; dim < N; ++dim)
        m_layer_tiles[dim] = LayerTilesType(init);
    }
    template <std::integral SizeType>
    Tiles(std::array<SizeType, N> sizes) {
      for (auto dim = 0u; dim < N; ++dim)
        m_layer_tiles[dim] = LayerTilesType(sizes[dim]);
    }

    const auto& operator[](std::size_t idx) const { return m_layer_tiles[idx]; }
    auto& operator[](std::size_t idx) { return m_layer_tiles[idx]; }

  private:
    std::array<LayerTiles, N> m_layer_tiles;
  };

}  // namespace ticl
