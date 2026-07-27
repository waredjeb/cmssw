// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <cstdint>
#include <numeric>
#include <span>
#include <vector>

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

    // Compute the maximum layer index over all rechits into d_max[0].
    struct MaxLayerKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                    int* d_max,
                                    const uint32_t size) const {
        for (auto i : uniform_elements(acc, size)) {
          alpaka::atomicAdd(acc, d_max, 0);  // no-op to keep acc used on all backends
          alpaka::atomicMax(acc, d_max, inputs[i].layer());
        }
      }
    };

    // Histogram: count rechits per layer.
    struct HistogramKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                    int* counts,
                                    const uint32_t size) const {
        for (auto i : uniform_elements(acc, size)) {
          alpaka::atomicAdd(acc, &counts[inputs[i].layer()], 1);
        }
      }
    };

    // Counting-sort scatter: build the permutation sorted-position -> original
    // rechit index using per-layer cursors (initialised to the layer offsets).
    struct PermKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                    int* cursors,
                                    int* perm,
                                    const uint32_t size) const {
        for (auto i : uniform_elements(acc, size)) {
          const int pos = alpaka::atomicAdd(acc, &cursors[inputs[i].layer()], 1);
          perm[pos] = static_cast<int>(i);
        }
      }
    };

    // Gather input columns into layer-contiguous buffers, in sorted order.
    struct GatherKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                    const int* perm,
                                    float* dim1,
                                    float* dim2,
                                    float* energy,
                                    float* sigma,
                                    const uint32_t size) const {
        for (auto s : uniform_elements(acc, size)) {
          const int orig = perm[s];
          dim1[s] = inputs[orig].dim1();
          dim2[s] = inputs[orig].dim2();
          energy[s] = inputs[orig].energy();
          sigma[s] = inputs[orig].sigmaNoise();
        }
      }
    };

    // Scatter cluster indices back to the original rechit order and reset seeds.
    struct ScatterKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int* perm,
                                    const int* clusterIndexSorted,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs,
                                    const uint32_t size) const {
        for (auto s : uniform_elements(acc, size)) {
          const int orig = perm[s];
          const int ci = clusterIndexSorted[s];
          outputs[orig].clusterIndex() = (ci < 0) ? kInvalidCluster : ci;
          outputs[orig].isSeed() = 0;
        }
      }
    };

    // Flag the seeds (seed indices are in sorted order -> map through perm).
    struct SeedKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int* perm,
                                    const int32_t* seeds,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs,
                                    const uint32_t nseeds) const {
        for (auto k : uniform_elements(acc, nseeds)) {
          const int sorted = seeds[k];
          outputs[perm[sorted]].isSeed() = 1;
        }
      }
    };

    // Reduce the maximum global cluster index into d_max[0].
    struct MaxClusterKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int* clusterIndexSorted,
                                    int* d_max,
                                    const uint32_t size) const {
        for (auto s : uniform_elements(acc, size)) {
          alpaka::atomicMax(acc, d_max, clusterIndexSorted[s]);
        }
      }
    };

    // Write the number of clusters (max global cluster index + 1) into the scalar.
    struct SetScalarKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int* d_max,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs) const {
        if (once_per_grid(acc)) {
          outputs.numberOfClustersScalar() = (d_max[0] < 0) ? 0u : static_cast<unsigned int>(d_max[0] + 1);
        }
      }
    };

  }  // namespace

  void HGCalCLUEsteringAlgoWrapper::run(Queue& queue,
                                        const unsigned int size,
                                        const float dc,
                                        const float kappa,
                                        const float outlierDeltaFactor,
                                        const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                        HGCalSoARecHitsExtraDeviceCollection::View outputs) const {
    // Nothing to do for an empty event: just publish zero clusters.
    if (size == 0) {
      auto nClusters = make_device_view<unsigned int>(queue, outputs.numberOfClustersScalar());
      alpaka::memset(queue, nClusters, 0x0);
      return;
    }

    const uint32_t items = 256;
    const uint32_t groups = divide_up_by(size, items);
    const auto workDiv = make_workdiv<Acc1D>(groups, items);

    // 1) Determine the number of layers present: L = max(layer) + 1.
    auto d_maxLayer = make_device_buffer<int[]>(queue, 1u);
    alpaka::memset(queue, d_maxLayer, 0x0);
    alpaka::exec<Acc1D>(queue, workDiv, MaxLayerKernel{}, inputs, d_maxLayer.data(), size);
    auto h_maxLayer = make_host_buffer<int[]>(queue, 1u);
    alpaka::memcpy(queue, h_maxLayer, d_maxLayer);
    alpaka::wait(queue);
    const int nLayers = h_maxLayer[0] + 1;

    // 2) Histogram the rechits per layer.
    auto d_counts = make_device_buffer<int[]>(queue, nLayers);
    alpaka::memset(queue, d_counts, 0x0);
    alpaka::exec<Acc1D>(queue, workDiv, HistogramKernel{}, inputs, d_counts.data(), size);

    // Copy counts to host and build the exclusive-scan offsets (= per-layer
    // starting positions) and the per-layer event sizes for the batched run.
    auto h_counts = make_host_buffer<int[]>(queue, nLayers);
    alpaka::memcpy(queue, h_counts, d_counts);
    alpaka::wait(queue);

    // event_sizes holds only the NON-EMPTY layers. CLUEstering's batched
    // clustering launches per-batch device work, and a zero-size batch triggers a
    // 0-block kernel launch on CUDA (cudaErrorInvalidValue; harmless no-op on the
    // serial backend). Empty layers contribute no points, so dropping them leaves
    // the contiguous per-layer point layout and the global cluster numbering
    // unchanged. h_offsets is still kept per layer for the sort/gather kernels.
    std::vector<uint32_t> event_sizes;
    event_sizes.reserve(nLayers);
    auto h_offsets = make_host_buffer<int[]>(queue, nLayers);
    int running = 0;
    for (int l = 0; l < nLayers; ++l) {
      h_offsets[l] = running;
      if (h_counts[l] > 0)
        event_sizes.push_back(static_cast<uint32_t>(h_counts[l]));
      running += h_counts[l];
    }

    // 3) Counting-sort: build permutation sorted-position -> original index.
    auto d_cursors = make_device_buffer<int[]>(queue, nLayers);
    alpaka::memcpy(queue, d_cursors, h_offsets);
    auto d_perm = make_device_buffer<int[]>(queue, size);
    alpaka::exec<Acc1D>(queue, workDiv, PermKernel{}, inputs, d_cursors.data(), d_perm.data(), size);

    // 4) Gather the input columns into layer-contiguous device buffers.
    auto d_dim1 = make_device_buffer<float[]>(queue, size);
    auto d_dim2 = make_device_buffer<float[]>(queue, size);
    auto d_energy = make_device_buffer<float[]>(queue, size);
    auto d_sigma = make_device_buffer<float[]>(queue, size);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        GatherKernel{},
                        inputs,
                        d_perm.data(),
                        d_dim1.data(),
                        d_dim2.data(),
                        d_energy.data(),
                        d_sigma.data(),
                        size);

    // Output buffer for the per-point (sorted-order) global cluster indices.
    auto d_clusterIndexSorted = make_device_buffer<int[]>(queue, size);
    alpaka::memset(queue, d_clusterIndexSorted, kInvalidClusterByte);
    alpaka::wait(queue);

#ifdef HGCAL_CLUESTERING_DEVICE_ENABLED
    // 5) Build the CLUEstering device points from the sorted buffers.
    //    Variadic ctor: Ndim coordinate buffers + weight buffer + int output.
    clue::PointsDevice<2> d_points(
        queue, static_cast<int32_t>(size), d_dim1.data(), d_dim2.data(), d_energy.data(), d_clusterIndexSorted.data());
    d_points.set_density_uncertainty(std::span<float>(d_sigma.data(), size));

    // 6) Run the batched clustering: cluster indexes come out GLOBAL across
    //    layers, outliers are -1.
    // clue::Clusterer(density_radius, min_density, outlier_distance): the third
    // argument is an absolute distance. CLUE defines the outlier distance as
    // outlierDeltaFactor * dc (== the CPU algo's deltao, e.g. 2.0 * 1.3 = 2.6),
    // so scale it here rather than passing the bare factor.
    clue::Clusterer<2> algo(queue, dc, kappa, dc * outlierDeltaFactor);
    algo.make_clusters(queue, d_points, std::span<const uint32_t>(event_sizes));
    alpaka::wait(queue);
#else
    // HIP fallback: no device clustering (see note at the top of the file).
    // d_clusterIndexSorted stays kInvalidCluster (-1) from the memset above.
    (void)dc;
    (void)kappa;
    (void)outlierDeltaFactor;
    (void)d_dim1;
    (void)d_dim2;
    (void)d_energy;
    (void)d_sigma;
    (void)event_sizes;
#endif

    // 7) Scatter results back to the original rechit order.
    alpaka::exec<Acc1D>(
        queue, workDiv, ScatterKernel{}, d_perm.data(), d_clusterIndexSorted.data(), outputs, size);

#ifdef HGCAL_CLUESTERING_DEVICE_ENABLED
    // Flag the seeds (seed indices reference sorted-order positions).
    std::span<const int32_t> seeds = algo.getSeeds();
    const uint32_t nseeds = static_cast<uint32_t>(seeds.size());
    if (nseeds > 0) {
      const uint32_t seedGroups = divide_up_by(nseeds, items);
      const auto seedWorkDiv = make_workdiv<Acc1D>(seedGroups, items);
      alpaka::exec<Acc1D>(queue, seedWorkDiv, SeedKernel{}, d_perm.data(), seeds.data(), outputs, nseeds);
    }
#endif

    // 8) Number of clusters = max global cluster index + 1 (0 if none).
    auto d_maxCluster = make_device_buffer<int[]>(queue, 1u);
    alpaka::memset(queue, d_maxCluster, kInvalidClusterByte);  // -1
    alpaka::exec<Acc1D>(queue, workDiv, MaxClusterKernel{}, d_clusterIndexSorted.data(), d_maxCluster.data(), size);
    alpaka::exec<Acc1D>(queue, make_workdiv<Acc1D>(1u, 1u), SetScalarKernel{}, d_maxCluster.data(), outputs);
    alpaka::wait(queue);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
