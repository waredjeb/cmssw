// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "HGCalLayerClustersAlgoWrapper.h"

#include "CLUEstering/CLUEstering.hpp"

#include <cmath>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  // CLUEstering reports the seeds as one hit index per cluster, and keeps its
  // per-hit seed flag in a buffer of its own. Copy that information over to the
  // column the downstream modules read, and publish the number of clusters.
  class HGCalLayerClustersMarkSeedsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  const int32_t* seeds,
                                  const unsigned int numberOfClusters,
                                  HGCalSoARecHitsExtraDeviceCollection::View outputs) const {
      if (once_per_grid(acc)) {
        outputs.numberOfClustersScalar() = numberOfClusters;
      }
      for (auto cluster : uniform_elements(acc, numberOfClusters)) {
        outputs[seeds[cluster]].isSeed() = 1;
      }
    }
  };

  void HGCalLayerClustersAlgoWrapper::run(Queue& queue,
                                          const unsigned int size,
                                          const float dc,
                                          const float kappa,
                                          const float outlierDeltaFactor,
                                          const bool isScintillator,
                                          std::span<const uint32_t> batchItemSizes,
                                          const HGCalSoARecHitsDeviceCollection::ConstView inputs,
                                          HGCalSoARecHitsExtraDeviceCollection::View outputs) const {
    // Dont allow empty inputs
    if (size == 0) {
      return;
    }

    clue::ConstPointsDevice<2, float> points(queue,
                                             size,
                                             inputs.dim1().data(),
                                             inputs.dim2().data(),
                                             inputs.energy().data(),
                                             outputs.clusterIndex().data());

    // kappa * sigmaNoise for each hit.
    points.set_density_uncertainty({inputs.sigmaNoise().data(), size});

    // Following Params of CLUE
    clue::Clusterer<2> clusterer(queue, dc, kappa, outlierDeltaFactor * dc, dc);

    // One batch item per non-empty layer. CLUEstering has no notion of layers
    // and takes the layer of a hit from its position in the collection, which is
    // why the RecHits SoA is filled grouped by layer.
    if (isScintillator) {
      // Scintillator cells cluster in (eta, phi): phi is periodic, with the
      // coordinate stored in [0, 2pi) as required by the periodic metric.
      constexpr float kTwoPi = 2.f * static_cast<float>(M_PI);
      clusterer.setWrappedCoordinates(0, 1);
      clusterer.make_clusters(queue,
                              points,
                              batchItemSizes,
                              clue::PeriodicEuclideanMetric<2, float>{0.f, kTwoPi},
                              clue::FlatKernel<float>{0.5f});
    } else {
      clusterer.make_clusters(queue,
                              points,
                              batchItemSizes,
                              clue::EuclideanMetric<2, float>{},
                              clue::FlatKernel<float>{0.5f});
    }

    // Same number of seeds as clusters
    const auto seeds = clusterer.getSeeds();
    const auto numberOfClusters = static_cast<unsigned int>(points.n_clusters());
    if (numberOfClusters > 0) {
      const uint32_t items = 64;
      const uint32_t groups = divide_up_by(numberOfClusters, items);
      auto workDiv = make_workdiv<Acc1D>(groups, items);
      alpaka::exec<Acc1D>(
          queue, workDiv, HGCalLayerClustersMarkSeedsKernel{}, seeds.data(), numberOfClusters, outputs);
    }

    //wait for the kernel above to be done reading them.
    alpaka::wait(queue);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
