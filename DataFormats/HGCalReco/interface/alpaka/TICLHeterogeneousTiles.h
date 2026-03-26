
#pragma once

#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/FillAssociator.h"
#include "DataFormats/HGCalReco/interface/LayerTileConcept.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include <alpaka/alpaka.hpp>
#include <xtd/xtd.h>
#include <concepts>
#include <cstdint>
#include <span>

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  namespace concepts = ::ticl::concepts;

  template <concepts::LayerTile T>
  class LayerTilesView {
  private:
    ::ticl::AssociationMapView<> associations;

  public:
    constexpr LayerTilesView() = default;
    constexpr LayerTilesView(::ticl::AssociationMapView<> assoc_view) : associations{assoc_view} {}

    ALPAKA_FN_ACC auto phiBin(float phi) const {
      const auto normPhi = normalizedPhi(phi);
      const auto r = T::nPhiBins * M_1_PI * 0.5f;
      const auto phiBin = (normPhi + M_PI) * r;

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

    // ALPAKA_FN_ACC auto operator[](std::size_t idx) const { return associations[idx]; }
    ALPAKA_FN_ACC auto operator[](std::size_t idx) { return associations[idx]; }

    ALPAKA_FN_ACC auto searchBox(float etaMin, float etaMax, float phiMin, float phiMax) const {
      // if (!std::is_same_v<T, ticl::TileConstantsBarrel>) {
      //   if (etaMin * etaMax < 0) {
      //     return std::array<int, 4>({{0, 0, 0, 0}});
      //   }
      // }
      if (etaMax - etaMin < 0) {
        return std::array<int, 4>({{0, 0, 0, 0}});
      }
      const auto etaBinMin = etaBin(etaMin);
      const auto etaBinMax = etaBin(etaMax);
      const auto phiBinMin = phiBin(phiMin);
      const auto phiBinMax = phiBin(phiMax);
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
    template <concepts::LayerTile T>
    ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                  LayerTilesView<T> tiles,
                                  std::span<float> etas,
                                  std::span<float> phis,
                                  uint32_t* associations) const {
      for (auto idx : alpaka::uniformElements(acc, etas.size())) {
        associations[idx] = tiles.globalBin(etas[idx], phis[idx]);
      }
    }
  };

  // TODO: check consistensy of integer signedness
  template <concepts::LayerTile T>
  class LayerTiles {
  public:
    using View = LayerTilesView<T>;

    LayerTiles(Queue& queue, std::integral auto nvalues, std::integral auto nkeys)
        : m_associations(queue, nvalues, nkeys) {}
    LayerTiles(std::integral auto nvalues, std::integral auto nkeys)
      requires std::same_as<Device, alpaka::DevCpu>
        : m_associations(cms::alpakatools::host(), nvalues, nkeys) {}

    ALPAKA_FN_HOST void fill(Queue& queue,
                             std::span<float> etas,
                             std::span<float> phis,
                             std::span<uint32_t> layer_cluster_ids) {
      auto associations = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, etas.size());

      const auto blocksize = 1024u;
      const auto gridsize = cms::alpakatools::divide_up_by(etas.size(), blocksize);
      auto work_division = cms::alpakatools::make_workdiv(gridsize, blocksize);
      alpaka::exec<Acc1D>(queue, work_division, KernelTilesAssociations{}, etas, phis, associations.data());

      associator::fill<Acc1D>(
          queue, m_associations.view(), layer_cluster_ids, std::span<uint32_t>(associations.data(), etas.size()));
    }

    ALPAKA_FN_HOST auto view() const { return View(m_associations.view()); }

    // TODO: do we want to be able to clear them?
    void clear(Queue& queue);

  private:
    PortableCollection<Device, ::ticl::AssociationMap<>> m_associations;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl

namespace ticl {

  template <concepts::LayerTile T>
  using LayerTilesHost = alpaka_serial_sync::LayerTiles<T>;

}  // namespace ticl
