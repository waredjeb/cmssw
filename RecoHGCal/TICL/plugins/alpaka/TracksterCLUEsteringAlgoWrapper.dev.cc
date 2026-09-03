#include <algorithm>
#include <cstdint>
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

    // Cylinder distance: max(transverse Euclidean / mean angular width ratio, |delta layer id|),
    // so the unit shape is a disc times an interval along the layer axis.
    //
    // The per-point transverse angular ratio is owned by the metric, not CLUEstering, keeping the
    // tile search +- dc per axis, not scaled
    //
    // ratios <= 1, or CLUEstering box search incorrect
    struct CylinderMetric {
      using value_type = float;

      const float* ratio;

      template <typename TView>
      ALPAKA_FN_HOST_ACC float operator()(const TView& points, std::size_t i, std::size_t j) const {
        const auto pi = points[static_cast<int>(i)];
        const auto pj = points[static_cast<int>(j)];
        const float dx = pi[0] - pj[0];
        const float dy = pi[1] - pj[1];
        const float dz = pi[2] - pj[2];
        const float transverse = clue::math::sqrt(dx * dx + dy * dy) * 2.f / (ratio[i] + ratio[j]);
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
    struct GatherPointsKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::reco::CaloClusterSoAConstView input,
                                    const int32_t* flags,
                                    const int32_t* scan,
                                    uint32_t offset,
                                    uint32_t size,
                                    float invSigmaRef,
                                    float sigmaRatio,
                                    float invLayerScale,
                                    float endcapGap,
                                    float rhocEtaExponent,
                                    float rhocPivotRadius,
                                    float* x,
                                    float* y,
                                    float* z,
                                    float* energy,
                                    float* sigma,
                                    float* rhocScale,
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
          // Gnomonic direction coordinates in units of the reference (largest) sigmaT, and the
          // layer index in units of layerScale. The SoA layer already runs 0..N-1 for z<0 and
          // N..2N-1 for z>0; the extra gap for z>0 puts the two endcaps out of reach of each
          // other in the third coordinate while keeping the z range, and so the tile size
          // CLUEstering derives from it, comparable to the search radius.
          const auto invAbsZ = 1.f / clue::math::fabs(pz);
          x[slot] = px * invAbsZ * invSigmaRef;
          y[slot] = py * invAbsZ * invSigmaRef;
          z[slot] = static_cast<float>(input.position()[i].layer()) * invLayerScale + ((pz > 0.f) ? endcapGap : 0.f);
          // Per-point transverse scale, as a ratio to the reference sigma (<= 1): the metric
          // divides pair distances by the mean ratio, giving each sub-detector its own sigmaT
          // without displacing same-direction points across sub-detector boundaries.
          sigma[slot] = sigmaRatio;
          // Per-point multiplier on rhoc: (r0 / r)^alpha with r the gnomonic radius, so the
          // seed threshold rises towards high eta where the pile-up density under a fixed
          // angular disc is largest. alpha = 0 keeps the threshold flat.
          const auto r = clue::math::sqrt(px * px + py * py) * invAbsZ;
          rhocScale[slot] = clue::math::pow(rhocPivotRadius / r, rhocEtaExponent);
          energy[slot] = input.energy()[i].energy();
          sourceIndex[slot] = g;
          // The seed DetId is used to break ties deterministically inside CLUE.
          tags[slot] = input.indexes()[i].seedID().rawId();
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
    auto sigma = make_device_buffer<float[]>(queue, nSurviving);
    auto rhocScale = make_device_buffer<float[]>(queue, nSurviving);
    auto sourceIndex = make_device_buffer<uint32_t[]>(queue, nSurviving);
    auto tags = make_device_buffer<uint32_t[]>(queue, nSurviving);
    auto clusterIndex = make_device_buffer<int32_t[]>(queue, nSurviving);

    // Coordinates are stored in units of the largest sigmaT, so the per-point ratios stay <= 1
    // and the metric can only inflate distances with respect to the stored coordinates: CLUE's
    // coordinate search box (+-dc per axis) then stays exact.
    float sigmaRef = 0.f;
    for (auto s : parameters.sigmaT)
      sigmaRef = std::max(sigmaRef, s);

    // make sure first layer of next end cap is furrtehr than the max reach of the last layer
    // of previous, either d_c or outlierDistance
    const float endcapGap = 2.f * std::max(parameters.dc, parameters.outlierDistance);

    for (size_t d = 0; d < inputs.size(); ++d) {
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
                          parameters.sigmaT[d] / sigmaRef,
                          1.f / parameters.layerScale,
                          endcapGap,
                          parameters.rhocEtaExponent,
                          parameters.rhocPivotRadius,
                          x.data(),
                          y.data(),
                          z.data(),
                          energy.data(),
                          sigma.data(),
                          rhocScale.data(),
                          sourceIndex.data(),
                          tags.data());
    }

    clue::PointsDevice<3, float> points(
        queue, static_cast<int32_t>(nSurviving), x.data(), y.data(), z.data(), energy.data(), clusterIndex.data());
    points.set_tags(std::span<const uint32_t>(tags.data(), nSurviving));
    // Per-point multiplier on rhoc, applied by CLUE in the seed condition and in the
    // seeding-vs-outlier distance switch
    points.set_density_uncertainty(std::span<const float>(rhocScale.data(), nSurviving));

    clue::Clusterer<3> algo(
        queue, parameters.dc, parameters.rhoc, parameters.outlierDistance, parameters.seedingDistance);
    algo.make_clusters(queue, points, CylinderMetric{sigma.data()});

    {
      const auto workDiv = make_workdiv<Acc1D>(divide_up_by(nSurviving, items), items);
      alpaka::exec<Acc1D>(
          queue, workDiv, ScatterAssignmentKernel{}, clusterIndex.data(), sourceIndex.data(), assignment, nSurviving);
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
