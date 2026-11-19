#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include <alpaka/alpaka.hpp>

#include "DataFormats/DetId/interface/DetId.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "HGCalCLUEsteringAlgoWrapper.h"

#include "CLUEstering/core/Clusterer.hpp"
#include "CLUEstering/core/ConvolutionalKernel.hpp"
#include "CLUEstering/core/DistanceMetrics.hpp"
#include "CLUEstering/data_structures/PointsDevice.hpp"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace {

    // Per seed: its rank in rechit order is the index of its layer cluster. seedEnd is the inclusive scan of isSeed.
    struct KernelRankSeeds {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    const int32_t* isSeed,
                                    const int32_t* seedEnd,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
                                    int32_t* rankedIndex,
                                    uint32_t nRecHits) const {
        for (auto i : uniform_elements(acc, nRecHits)) {
          if (isSeed[i]) {
            rankedIndex[recHitsClustering[i].clusterIndex()] = seedEnd[i] - 1;
          }
        }
      }
    };

    // Per rechit: its layer cluster as the rank of the seed, its seed flag, and the number of layer clusters
    struct KernelMarkSeeds {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    const int32_t* isSeed,
                                    const int32_t* rankedIndex,
                                    HGCalSoARecHitsExtraDeviceCollection::View recHitsClustering,
                                    uint32_t nRecHits) const {
        for (auto i : uniform_elements(acc, nRecHits)) {
          const int cluster = recHitsClustering[i].clusterIndex();
          if (cluster >= 0) {
            recHitsClustering[i].clusterIndex() = rankedIndex[cluster];
          }
          recHitsClustering[i].isSeed() = isSeed[i];
          if (isSeed[i]) {
            alpaka::atomicAdd(acc, &recHitsClustering.numberOfClustersScalar(), 1u);
          }
        }
      }
    };

    // Per rechit: its layer if the layer clustering uses it, -1 otherwise
    struct KernelSelectRecHits {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    HGCalCLUEsteringAlgoWrapper::Selection selection,
                                    int32_t* layerKey,
                                    uint32_t nRecHits) const {
        for (auto i : uniform_elements(acc, nRecHits)) {
          const auto recHit = recHits[i];
          // DetId::det(), which is not available in device code
          const uint32_t det = (recHit.detid() >> DetId::kDetOffset) & DetId::kDetMask;
          const float energy = recHit.energy();
          const float sigmaNoise = recHit.sigmaNoise();
          const bool selected = ((selection.detMask >> det) & 1u) and recHit.layer() >= 0 and
                                std::isfinite(recHit.dim1()) and std::isfinite(recHit.dim2()) and
                                std::isfinite(energy) and std::isfinite(sigmaNoise) and sigmaNoise >= 0.f and
                                energy > selection.minEnergy and energy > selection.ecut * sigmaNoise;
          layerKey[i] = selected ? recHit.layer() : -1;
        }
      }
    };

    // Per selected rechit, in layer order: copies the rechit from the full collection
    struct KernelGatherRecHits {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    const int32_t* order,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    HGCalSoARecHitsDeviceCollection::View selected,
                                    uint32_t nSelected) const {
        for (auto j : uniform_elements(acc, nSelected)) {
          const auto in = recHits[order[j]];
          auto out = selected[j];
          out.dim1() = in.dim1();
          out.dim2() = in.dim2();
          out.dim3() = in.dim3();
          out.layer() = in.layer();
          out.energy() = in.energy();
          out.mipEnergy() = in.mipEnergy();
          out.sigmaNoise() = in.sigmaNoise();
          out.recHitIndex() = in.recHitIndex();
          out.flags() = in.flags();
          out.detid() = in.detid();
          out.time() = in.time();
          out.timeError() = in.timeError();
        }
      }
    };

    // Per selected rechit: moves its layer cluster and seed flag back to its row of the full collection
    struct KernelScatterClustering {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    const int32_t* order,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView selectedClustering,
                                    HGCalSoARecHitsExtraDeviceCollection::View recHitsClustering,
                                    uint32_t nSelected) const {
        if (once_per_grid(acc)) {
          recHitsClustering.numberOfClustersScalar() = selectedClustering.numberOfClustersScalar();
        }
        for (auto j : uniform_elements(acc, nSelected)) {
          recHitsClustering[order[j]].clusterIndex() = selectedClustering[j].clusterIndex();
          recHitsClustering[order[j]].isSeed() = selectedClustering[j].isSeed();
        }
      }
    };

  }  // namespace

  struct HGCalCLUEsteringAlgoWrapper::Clusterer {
    clue::Clusterer<2> algo;
  };

  HGCalCLUEsteringAlgoWrapper::HGCalCLUEsteringAlgoWrapper(float dc,
                                                           float kappa,
                                                           float outlierDeltaFactor,
                                                           bool isScintillator)
      : clusterer_{std::make_unique<Clusterer>(Clusterer{clue::Clusterer<2>(dc, kappa, dc * outlierDeltaFactor)})},
        isScintillator_{isScintillator} {
    if (isScintillator_) {
      clusterer_->algo.setWrappedCoordinates(0, 1);
    }
  }

  HGCalCLUEsteringAlgoWrapper::~HGCalCLUEsteringAlgoWrapper() = default;

  void HGCalCLUEsteringAlgoWrapper::setExtremes(std::pair<float, float> dim1, std::pair<float, float> dim2) {
    const std::pair<float, float> extremes[2] = {dim1, dim2};
    clusterer_->algo.setExtremes(extremes);
  }

  void HGCalCLUEsteringAlgoWrapper::run(Queue& queue,
                                        uint32_t nRecHits,
                                        std::span<const uint32_t> batchItemSizes,
                                        HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                        HGCalSoARecHitsExtraDeviceCollection::View recHitsClustering) {
    auto nClusters = make_device_view<unsigned int>(queue, recHitsClustering.numberOfClustersScalar());
    alpaka::memset(queue, nClusters, 0x0);
    if (nRecHits == 0) {
      return;
    }

    clue::ConstPointsDevice<2, float> d_points(queue,
                                               static_cast<int32_t>(nRecHits),
                                               recHits.dim1().data(),
                                               recHits.dim2().data(),
                                               recHits.energy().data(),
                                               recHitsClustering.clusterIndex().data());
    d_points.set_density_uncertainty(recHits.sigmaNoise());
    d_points.set_tags(recHits.detid());

    if (isScintillator_) {
      constexpr float kTwoPi = 2.f * static_cast<float>(M_PI);
      clusterer_->algo.make_clusters(queue,
                                     d_points,
                                     batchItemSizes,
                                     clue::PeriodicEuclideanMetric<2, float>{0.f, kTwoPi},
                                     clue::FlatKernel<float>{0.5f});
    } else {
      clusterer_->algo.make_clusters(queue, d_points, batchItemSizes);
    }

    // The layer clusters are numbered by their seeds in rechit order, so the result does not depend on the order in
    // which CLUEstering finds the seeds. d_points releases is_seed only after the queued kernels that read it.
    const int32_t* isSeed = d_points.view().is_seed().data();
    auto seedEnd = make_device_buffer<int32_t[]>(queue, nRecHits);
    iterativePrefixScan<Acc1D>(queue, isSeed, seedEnd.data(), nRecHits);
    // there are at most as many layer clusters as rechits
    auto rankedIndex = make_device_buffer<int32_t[]>(queue, nRecHits);
    constexpr uint32_t blockSize = 256;
    const auto workDivRecHits = make_workdiv<Acc1D>(divide_up_by(nRecHits, blockSize), blockSize);
    alpaka::exec<Acc1D>(queue,
                        workDivRecHits,
                        KernelRankSeeds{},
                        isSeed,
                        seedEnd.data(),
                        recHitsClustering,
                        rankedIndex.data(),
                        nRecHits);
    alpaka::exec<Acc1D>(
        queue, workDivRecHits, KernelMarkSeeds{}, isSeed, rankedIndex.data(), recHitsClustering, nRecHits);
  }

  void HGCalCLUEsteringAlgoWrapper::runOnSelection(Queue& queue,
                                                   uint32_t nRecHits,
                                                   Selection const& selection,
                                                   HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                                   HGCalSoARecHitsExtraDeviceCollection::View recHitsClustering) {
    // the rechits that are not selected have no layer cluster
    auto nClusters = make_device_view<unsigned int>(queue, recHitsClustering.numberOfClustersScalar());
    alpaka::memset(queue, nClusters, 0x0);
    if (nRecHits == 0) {
      return;
    }
    auto clusterIndexView = make_device_view<int>(queue, recHitsClustering.clusterIndex().data(), nRecHits);
    alpaka::memset(queue, clusterIndexView, 0xff);
    auto isSeedView = make_device_view<uint8_t>(queue, recHitsClustering.isSeed().data(), nRecHits);
    alpaka::memset(queue, isSeedView, 0x0);

    constexpr uint32_t blockSize = 256;
    auto layerKey = make_device_buffer<int32_t[]>(queue, nRecHits);
    alpaka::exec<Acc1D>(queue,
                        make_workdiv<Acc1D>(divide_up_by(nRecHits, blockSize), blockSize),
                        KernelSelectRecHits{},
                        recHits,
                        selection,
                        layerKey.data(),
                        nRecHits);
    auto h_layerKey = make_host_buffer<int32_t[]>(queue, nRecHits);
    alpaka::memcpy(queue, h_layerKey, layerKey);
    // CLUEstering takes the batch sizes on the host
    alpaka::wait(queue);

    // Counting sort by layer, stable in rechit order: the selection and the batches do not depend on the backend.
    // offsets[layer + 1] first counts the selected rechits of each layer, then, after the running sum, offsets[layer]
    // is the position of the first rechit of the layer in the selection.
    int32_t maxLayer = -1;
    for (uint32_t i = 0; i < nRecHits; ++i) {
      maxLayer = std::max(maxLayer, h_layerKey[i]);
    }
    if (maxLayer < 0) {
      return;
    }
    std::vector<uint32_t> offsets(maxLayer + 2, 0);
    for (uint32_t i = 0; i < nRecHits; ++i) {
      if (h_layerKey[i] >= 0) {
        ++offsets[h_layerKey[i] + 1];
      }
    }
    std::vector<uint32_t> batches;
    for (int32_t layer = 0; layer <= maxLayer; ++layer) {
      if (offsets[layer + 1] > 0) {
        batches.push_back(offsets[layer + 1]);
      }
      offsets[layer + 1] += offsets[layer];
    }
    const uint32_t nSelected = offsets[maxLayer + 1];
    // order[j]: the rechit of the full collection at position j of the selection
    auto h_order = make_host_buffer<int32_t[]>(queue, nSelected);
    for (uint32_t i = 0; i < nRecHits; ++i) {
      if (h_layerKey[i] >= 0) {
        h_order[offsets[h_layerKey[i]]++] = static_cast<int32_t>(i);
      }
    }
    auto order = make_device_buffer<int32_t[]>(queue, nSelected);
    alpaka::memcpy(queue, order, h_order);

    // the selected rechits, contiguous by layer as CLUEstering needs them
    HGCalSoARecHitsDeviceCollection selected(queue, nSelected);
    const auto workDivSelected = make_workdiv<Acc1D>(divide_up_by(nSelected, blockSize), blockSize);
    alpaka::exec<Acc1D>(
        queue, workDivSelected, KernelGatherRecHits{}, order.data(), recHits, selected.view(), nSelected);

    HGCalSoARecHitsExtraDeviceCollection selectedClustering(queue, nSelected);
    run(queue, nSelected, std::span<const uint32_t>(batches), selected.const_view(), selectedClustering.view());

    alpaka::exec<Acc1D>(queue,
                        workDivSelected,
                        KernelScatterClustering{},
                        order.data(),
                        selectedClustering.const_view(),
                        recHitsClustering,
                        nSelected);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
