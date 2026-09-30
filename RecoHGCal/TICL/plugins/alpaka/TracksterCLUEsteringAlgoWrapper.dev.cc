#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "TracksterCLUEsteringAlgoWrapper.h"

#include "CLUEstering/core/Clusterer.hpp"
#include "CLUEstering/data_structures/PointsDevice.hpp"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace {

    // Cylinder distance: max(transverse Euclidean / mean relative sigmaT of the pair, |delta layer id|),
    // so the unit shape is a disc times an interval along the layer axis.
    //
    // The per-point relative sigmaT is owned by the metric, not CLUEstering, keeping the
    // tile search +- dc per axis, not scaled
    //
    // relativeSigmaT >= 1, so risk of shrinking distances, so the box is widened per point through CLUEsterings
    // sigma.
    struct CylinderMetric {
      using value_type = float;

      // per point: sigmaT of its detector / smallest sigmaT
      const float* relativeSigmaT;

      template <typename TView>
      ALPAKA_FN_HOST_ACC float operator()(const TView& points, std::size_t i, std::size_t j) const {
        const auto pi = points[static_cast<int>(i)];
        const auto pj = points[static_cast<int>(j)];
        const float dx = pi[0] - pj[0];
        const float dy = pi[1] - pj[1];
        const float dz = pi[2] - pj[2];
        const float transverse = clue::math::sqrt(dx * dx + dy * dy) * 2.f / (relativeSigmaT[i] + relativeSigmaT[j]);
        return clue::math::max(transverse, clue::math::fabs(dz));
      }
    };

    // flags[i] = 1 for the layer clusters that survive the iteration mask, 0 otherwise.
    struct MaskToFlagsKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc, const float* mask, int32_t* flags, uint32_t size) const {
        for (auto i : uniform_elements(acc, size)) {
          flags[i] = (mask[i] > 0.f) ? 1 : 0;
        }
      }
    };

    // Copies the surviving layer clusters of one detector into the compacted point buffers.
    // `scan` is the inclusive prefix sum of the flags over the whole merged collection, so
    // scan[g] - 1 is the slot of the layer cluster at merged index g.
    // `layerRange` is the [min, max] layer of points in the detector
    // as {min-, max-, min+, max+}
    struct GatherPointsKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::reco::CaloClusterSoAConstView input,
                                    const int32_t* flags,
                                    const int32_t* scan,
                                    uint32_t offset,
                                    uint32_t size,
                                    float invSigmaRef,
                                    float detectorRelativeSigmaT,
                                    float invLayerScale,
                                    float endcapGap,
                                    float rhocEtaExponent,
                                    float rhocPivotRadius,
                                    float* x,
                                    float* y,
                                    float* z,
                                    float* energy,
                                    float* relativeSigmaT,
                                    float* rhocScale,
                                    int32_t* layer,
                                    uint8_t* side,
                                    int32_t* layerRange,
                                    uint32_t* sourceIndex,
                                    uint32_t* tags) const {
        for (auto i : uniform_elements(acc, size)) {
          const auto g = offset + i;
          if (flags[g] == 0)
            continue;
          const auto slot = static_cast<uint32_t>(scan[g]) - 1u;
          const auto px = input.position()[i].x();
          const auto py = input.position()[i].y();
          const auto pz = input.position()[i].z();
          const auto layerId = static_cast<int32_t>(input.position()[i].layer());
          const auto sideId = static_cast<uint8_t>(pz > 0.f);

          // Gnomonic direction coordinates in units of the smallest sigmaT, with z as the layer index
          const auto invAbsZ = 1.f / clue::math::fabs(pz);
          x[slot] = px * invAbsZ * invSigmaRef;
          y[slot] = py * invAbsZ * invSigmaRef;
          z[slot] = static_cast<float>(layerId) * invLayerScale + (sideId ? endcapGap : 0.f);

          // same for every point of the detector, calculated on host
          relativeSigmaT[slot] = detectorRelativeSigmaT;
          layer[slot] = layerId;
          side[slot] = sideId;
          alpaka::atomicMin(acc, layerRange + 2 * sideId, layerId);
          alpaka::atomicMax(acc, layerRange + 2 * sideId + 1, layerId);

          // rhoc scales with (r0 / r) ^ alpha, with r the gnomonic radius
          // r = 0 would give infinity, but this shouldn't be possible
          const auto r = clue::math::sqrt(px * px + py * py) * invAbsZ;
          rhocScale[slot] = clue::math::pow(rhocPivotRadius / r, rhocEtaExponent);
          energy[slot] = input.energy()[i].energy();
          sourceIndex[slot] = g;
          // The seed DetId is used to break ties deterministically inside CLUEstering.
          tags[slot] = input.indexes()[i].seedID().rawId();
        }
      }
    };

    // Per-point transverse search window. A neighbour j is within metric distance R of point i
    // only if its stored transverse distance is <= R * (relativeSigmaT_i + relativeSigmaT_j) / 2, so the widest
    // box point i needs is set by the largest relativeSigmaT among the detectors it can reach along the layer
    // axis. CLUEstering widens the box of dimension d to R * sigma_d * sqrt(2), hence the sigma.
    struct SearchWindowKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int32_t* layer,
                                    const uint8_t* side,
                                    const float* relativeSigmaT,
                                    const int32_t* layerRange,
                                    const float* detectorRelativeSigmaT,
                                    uint32_t nDetectors,
                                    int32_t layerReach,
                                    float* searchSigma,
                                    uint32_t size) const {
        // rounded up, so that sigma * sqrt(2) in CLUEstering is never below (relativeSigmaT_i + max) / 2
        constexpr float invSqrt2 = 0.70710681f;
        for (auto i : uniform_elements(acc, size)) {
          const auto low = layer[i] - layerReach;
          const auto high = layer[i] + layerReach;
          auto maxRelativeSigmaT = relativeSigmaT[i];
          for (uint32_t d = 0; d < nDetectors; ++d) {
            const auto* range = layerRange + 4 * d + 2 * side[i];
            if (range[0] <= high and range[1] >= low)
              maxRelativeSigmaT = clue::math::max(maxRelativeSigmaT, detectorRelativeSigmaT[d]);
          }
          searchSigma[i] = 0.5f * (relativeSigmaT[i] + maxRelativeSigmaT) * invSqrt2;
        }
      }
    };

    // Writes the trackster index found for each compacted point back to merged-LC indexing.
    struct ScatterAssignmentKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int32_t* clusterIndex,
                                    const uint32_t* sourceIndex,
                                    int32_t* assignment,
                                    uint32_t size) const {
        for (auto k : uniform_elements(acc, size)) {
          assignment[sourceIndex[k]] = clusterIndex[k];
        }
      }
    };

  }  // namespace

  void TracksterCLUEsteringAlgoWrapper::run(Queue& queue,
                                            std::span<const ::reco::CaloClusterSoAConstView> inputs,
                                            std::span<const uint32_t> offsets,
                                            const float* mask,
                                            uint32_t nTotal,
                                            uint32_t nSurviving,
                                            Parameters const& parameters,
                                            int32_t* assignment) const {
    // Masked layer clusters and CLUE outliers alike are left unassigned.
    auto assignmentView = make_device_view(queue, assignment, nTotal);
    alpaka::fill(queue, assignmentView, static_cast<int32_t>(-1));

    if (nTotal == 0 or nSurviving == 0)
      return;

    constexpr auto items = 256u;
    const auto nDetectors = static_cast<uint32_t>(inputs.size());

    auto flags = make_device_buffer<int32_t[]>(queue, nTotal);
    auto scan = make_device_buffer<int32_t[]>(queue, nTotal);
    {
      const auto workDiv = make_workdiv<Acc1D>(divide_up_by(nTotal, items), items);
      alpaka::exec<Acc1D>(queue, workDiv, MaskToFlagsKernel{}, mask, flags.data(), nTotal);
    }
    iterativePrefixScan<Acc1D>(flags.data(), scan.data(), nTotal, queue);

    auto x = make_device_buffer<float[]>(queue, nSurviving);
    auto y = make_device_buffer<float[]>(queue, nSurviving);
    auto z = make_device_buffer<float[]>(queue, nSurviving);
    auto energy = make_device_buffer<float[]>(queue, nSurviving);
    auto relativeSigmaT = make_device_buffer<float[]>(queue, nSurviving);
    auto rhocScale = make_device_buffer<float[]>(queue, nSurviving);
    auto layer = make_device_buffer<int32_t[]>(queue, nSurviving);
    auto side = make_device_buffer<uint8_t[]>(queue, nSurviving);
    auto searchSigma = make_device_buffer<float[]>(queue, nSurviving);
    auto sourceIndex = make_device_buffer<uint32_t[]>(queue, nSurviving);
    auto tags = make_device_buffer<uint32_t[]>(queue, nSurviving);
    auto clusterIndex = make_device_buffer<int32_t[]>(queue, nSurviving);

    // Coordinates are stored in units of the smallest sigmaT, so relativeSigmaT = sigmaT / sigmaRef >= 1:
    // the search box of a point is then at least +-dc in stored coordinates and is widened per
    // point by SearchWindowKernel to cover the coarser detectors it can reach.
    const float sigmaRef = std::ranges::min(parameters.sigmaT);

    // make sure first layer of next end cap is furrtehr than the max reach of the last layer
    // of previous, either d_c or outlierDistance
    const float endcapGap = 2.f * std::max(parameters.dc, parameters.outlierDistance);

    // Per detector: its relative sigmaT, and the [min, max] layer of its surviving points per side.
    auto detectorRelativeSigmaTHost = make_host_buffer<float[]>(queue, nDetectors);
    auto layerRangeHost = make_host_buffer<int32_t[]>(queue, 4 * nDetectors);
    for (uint32_t d = 0; d < nDetectors; ++d) {
      detectorRelativeSigmaTHost[d] = parameters.sigmaT[d] / sigmaRef;
      for (auto s = 0; s < 2; ++s) {
        layerRangeHost[4 * d + 2 * s] = std::numeric_limits<int32_t>::max();
        layerRangeHost[4 * d + 2 * s + 1] = std::numeric_limits<int32_t>::min();
      }
    }
    auto detectorRelativeSigmaT = make_device_buffer<float[]>(queue, nDetectors);
    auto layerRange = make_device_buffer<int32_t[]>(queue, 4 * nDetectors);
    alpaka::memcpy(queue, detectorRelativeSigmaT, detectorRelativeSigmaTHost);
    alpaka::memcpy(queue, layerRange, layerRangeHost);

    for (uint32_t d = 0; d < nDetectors; d++) {
      const auto size = static_cast<uint32_t>(inputs[d].metadata().size()[0]);
      if (size == 0)
        continue;
      const auto workDiv = make_workdiv<Acc1D>(divide_up_by(size, items), items);
      alpaka::exec<Acc1D>(queue,
                          workDiv,
                          GatherPointsKernel{},
                          inputs[d],
                          flags.data(),
                          scan.data(),
                          offsets[d],
                          size,
                          1.f / sigmaRef,
                          detectorRelativeSigmaTHost[d],
                          1.f / parameters.layerScale,
                          endcapGap,
                          parameters.rhocEtaExponent,
                          parameters.rhocPivotRadius,
                          x.data(),
                          y.data(),
                          z.data(),
                          energy.data(),
                          relativeSigmaT.data(),
                          rhocScale.data(),
                          layer.data(),
                          side.data(),
                          layerRange.data() + 4 * d,
                          sourceIndex.data(),
                          tags.data());
    }

    {
      // The largest search box along the layer axis is +-max(dc, outlierDistance), in layers.
      const auto layerReach =
          static_cast<int32_t>(std::ceil(std::max(parameters.dc, parameters.outlierDistance) * parameters.layerScale));
      const auto workDiv = make_workdiv<Acc1D>(divide_up_by(nSurviving, items), items);
      alpaka::exec<Acc1D>(queue,
                          workDiv,
                          SearchWindowKernel{},
                          layer.data(),
                          side.data(),
                          relativeSigmaT.data(),
                          layerRange.data(),
                          detectorRelativeSigmaT.data(),
                          nDetectors,
                          layerReach,
                          searchSigma.data(),
                          nSurviving);
    }

    clue::PointsDevice<3, float> points(
        queue, static_cast<int32_t>(nSurviving), x.data(), y.data(), z.data(), energy.data(), clusterIndex.data());
    points.set_tags(std::span<const uint32_t>(tags.data(), nSurviving));
    // Per-point multiplier on rhoc, applied by CLUE in the seed condition and in the
    // seeding-vs-outlier distance switch
    points.set_density_uncertainty(std::span<const float>(rhocScale.data(), nSurviving));
    // Per-point transverse search window; the layer axis keeps the plain +-R box.
    points.set_sigma(0, std::span<const float>(searchSigma.data(), nSurviving));
    points.set_sigma(1, std::span<const float>(searchSigma.data(), nSurviving));

    clue::Clusterer<3> algo(
        queue, parameters.dc, parameters.rhoc, parameters.outlierDistance, parameters.seedingDistance);
    algo.make_clusters(queue, points, CylinderMetric{relativeSigmaT.data()});

    {
      const auto workDiv = make_workdiv<Acc1D>(divide_up_by(nSurviving, items), items);
      alpaka::exec<Acc1D>(
          queue, workDiv, ScatterAssignmentKernel{}, clusterIndex.data(), sourceIndex.data(), assignment, nSurviving);
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
