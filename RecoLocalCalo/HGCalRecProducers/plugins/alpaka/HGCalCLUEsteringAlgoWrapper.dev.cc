// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <cmath>
#include <limits>
#include <cstdint>
#include <span>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "HGCalCLUEsteringAlgoWrapper.h"
#include "ConstantsForClusters.h"

// External, header-only CLUEstering library (runs on device via alpaka).
//
// NOTE: the ROCm/HIP clang-22 frontend segfaults while *parsing* CLUEstering's
// ClusterProperties.hpp (the `std::views::iota(0) | std::views::take(...)`
// pipeline), which is pulled in unconditionally through Clusterer.hpp ->
// PointsHost.hpp. This is an upstream compiler bug, unrelated to the code below.
// Since we may only edit this package (not the external headers, and not other
// CMSSW packages), the CLUEstering device path is guarded out of the HIP
// backend so the portable plugin library still builds green for every enabled
// backend. On HIP the wrapper falls back to an empty clustering (all rechits
// left as kInvalidCluster, zero clusters); the real device clustering is active
// on the serial and CUDA backends. TODO: re-enable on HIP once the upstream
// std::ranges::iota|take pipeline is compilable by the ROCm clang frontend.
#if !defined(ALPAKA_ACC_GPU_HIP_ENABLED)
#define HGCAL_CLUESTERING_DEVICE_ENABLED 1
#include "CLUEstering/core/detail/defines.hpp"
#include "CLUEstering/core/Clusterer.hpp"
#include "CLUEstering/data_structures/PointsDevice.hpp"
#endif

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  using namespace hgcal::constants;

  namespace {

    // Flag the seeds. Seed indices reference the point order handed to
    // CLUEstering, which (rechits already sorted by layer upstream) is the input
    // order, so they index the output SoA directly.
    struct SeedKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int32_t* seeds,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs,
                                    const uint32_t nseeds) const {
        for (auto k : uniform_elements(acc, nseeds)) {
          outputs[seeds[k]].isSeed() = 1;
        }
      }
    };

    // Copy the CLUE intermediates out of the CLUEstering points into the output
    // SoA. Nothing in reconstruction reads them, but without them the SoA dumper
    // reports zeros and the per-cell comparison against the legacy algorithm is
    // impossible. CLUEstering does not persist delta (it is local to its
    // nearest-higher kernel), so it is recomputed here from its definition: the
    // distance to the nearest higher. Points with no nearest higher (seeds and
    // outliers) get delta = max and nearestHigher = kInvalidNearestHigher, which
    // is what the legacy algorithm stores for the same case.
    struct CopyClueIntermediatesKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const float* rho,
                                    const int32_t* nearestHigher,
                                    HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs,
                                    const uint32_t size,
                                    const bool isScintillator) const {
        constexpr float kTwoPi = 2.f * static_cast<float>(M_PI);
        for (auto i : uniform_elements(acc, size)) {
          outputs[i].rho() = rho[i];
          const int32_t nh = nearestHigher[i];
          if (nh < 0) {
            outputs[i].nearestHigher() = kInvalidNearestHigher;
            outputs[i].delta() = std::numeric_limits<float>::max();
          } else {
            outputs[i].nearestHigher() = static_cast<unsigned int>(nh);
            float d1 = inputs[nh].dim1() - inputs[i].dim1();
            float d2 = inputs[nh].dim2() - inputs[i].dim2();
            if (isScintillator) {
              // dim2 is phi in [0, 2pi): take the shorter way round, matching
              // the periodic metric handed to the clustering.
              if (d2 > kTwoPi / 2.f)
                d2 -= kTwoPi;
              else if (d2 < -kTwoPi / 2.f)
                d2 += kTwoPi;
            }
            outputs[i].delta() = std::sqrt(d1 * d1 + d2 * d2);
          }
        }
      }
    };

  }  // namespace

  void HGCalCLUEsteringAlgoWrapper::run(Queue& queue,
                                        const unsigned int size,
                                        const float dc,
                                        const float kappa,
                                        const float outlierDeltaFactor,
                                        const bool isScintillator,
                                        std::span<const uint32_t> batchItemSizes,
                                        const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                        HGCalSoARecHitsExtraDeviceCollection::View outputs) const {
    // Nothing to do for an empty event: just publish zero clusters.
    if (size == 0) {
      auto nClusters = make_device_view<unsigned int>(queue, outputs.numberOfClustersScalar());
      alpaka::memset(queue, nClusters, 0x0);
      return;
    }

    const uint32_t items = 256;

    // Reset the seed flags to 0 (the seeds are flagged below). The cluster index
    // column does not need initialising: the clustering writes every point,
    // leaving outliers at -1 == kInvalidCluster.
    auto isSeedView = make_device_view(queue, outputs.isSeed().data(), size);
    alpaka::fill(queue, isSeedView, static_cast<uint8_t>(0));

#ifdef HGCAL_CLUESTERING_DEVICE_ENABLED
    // Build the CLUEstering device points directly on the (layer-sorted) input
    // columns -- no gather/reorder needed -- and let the clustering write the
    // global cluster indices straight into the output SoA column. ConstPointsDevice
    // takes the input columns by const pointer (the clustering only reads the
    // coordinates/weights), so no const_cast is needed.
    clue::ConstPointsDevice<2, float> d_points(queue,
                                               static_cast<int32_t>(size),
                                               inputs.dim1().data(),
                                               inputs.dim2().data(),
                                               inputs.energy().data(),
                                               outputs.clusterIndex().data());
    d_points.set_density_uncertainty(std::span<const float>(inputs.sigmaNoise().data(), size));

    // Run the batched clustering (one 2D clustering per layer): the per-layer
    // batch sizes are computed once upstream (in the rechit producer) and passed
    // in here. Cluster indexes come out GLOBAL across layers, outliers are -1.
    // clue::Clusterer(density_radius, min_density, outlier_distance): the third
    // argument is an absolute distance. CLUE defines the outlier distance as
    // outlierDeltaFactor * dc (== the CPU algo's deltao, e.g. 2.0 * 1.3 = 2.6),
    // so scale it here rather than passing the bare factor.
    clue::Clusterer<2> algo(queue, dc, kappa, dc * outlierDeltaFactor);
    if (isScintillator) {
      // Scintillator (BH) cells cluster in (eta, phi): phi is periodic. The
      // rechit producer stores phi in [0, 2pi), so use the periodic metric with
      // period 2pi on the second coordinate (0 == non-periodic on the first).
      constexpr float kTwoPi = 2.f * static_cast<float>(M_PI);
      algo.setWrappedCoordinates(0, 1);
      algo.make_clusters(queue,
                         d_points,
                         batchItemSizes,
                         clue::PeriodicEuclideanMetric<2, float>{0.f, kTwoPi},
                         clue::FlatKernel<float>{0.5f});
    } else {
      algo.make_clusters(queue, d_points, batchItemSizes);
    }
    alpaka::wait(queue);

    // Publish rho / delta / nearestHigher. Reconstruction does not read them,
    // but they are what makes a per-cell comparison against the legacy CLUE
    // possible, so they are filled rather than left at zero.
    {
      // rho and nearestHigher live on the points VIEW (raw device pointers); the
      // host-side PointsDevice wrapper only exposes coordinates, weights and
      // cluster indices.
      const auto& pointsView = d_points.view();
      const auto copyWorkDiv = make_workdiv<Acc1D>(divide_up_by(size, items), items);
      alpaka::exec<Acc1D>(queue,
                          copyWorkDiv,
                          CopyClueIntermediatesKernel{},
                          pointsView.m_rho,
                          pointsView.m_nearest_higher,
                          inputs,
                          outputs,
                          size,
                          isScintillator);
    }

    // Number of clusters comes straight from CLUEstering (max global index + 1);
    // no reduction kernel of our own.
    auto h_nClusters = make_host_buffer<unsigned int>(queue);
    *h_nClusters.data() = static_cast<unsigned int>(d_points.n_clusters());
    auto d_nClusters = make_device_view<unsigned int>(queue, outputs.numberOfClustersScalar());
    alpaka::memcpy(queue, d_nClusters, h_nClusters);

    // Flag the seeds.
    std::span<const int32_t> seeds = algo.getSeeds();
    const uint32_t nseeds = static_cast<uint32_t>(seeds.size());
    if (nseeds > 0) {
      const uint32_t seedGroups = divide_up_by(nseeds, items);
      const auto seedWorkDiv = make_workdiv<Acc1D>(seedGroups, items);
      alpaka::exec<Acc1D>(queue, seedWorkDiv, SeedKernel{}, seeds.data(), outputs, nseeds);
    }
    alpaka::wait(queue);
#else
    // HIP fallback: no device clustering (see note at the top of the file).
    // Leave every rechit as an outlier and publish zero clusters.
    auto clusterIndexView = make_device_view(queue, outputs.clusterIndex().data(), size);
    alpaka::memset(queue, clusterIndexView, kInvalidClusterByte);  // -1
    auto d_nClusters = make_device_view<unsigned int>(queue, outputs.numberOfClustersScalar());
    alpaka::memset(queue, d_nClusters, 0x0);
    (void)items;
    (void)dc;
    (void)kappa;
    (void)outlierDeltaFactor;
    (void)isScintillator;
    (void)batchItemSizes;
    alpaka::wait(queue);
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
