// Smoke test: verify the external CLUEstering (master) API surface used by the
// legacy HGCalCLUEAlgo port compiles/type-checks in this CMSSW/alpaka environment.
// It is not meant to produce meaningful clustering output.

#include "CLUEstering/CLUEstering.hpp"

#include <span>
#include <vector>
#include <ranges>
#include <iostream>

int main() {
  // A handful of 2D points with weights.
  std::vector<float> x{0.f, 0.1f, 5.f, 5.1f};
  std::vector<float> y{0.f, 0.1f, 5.f, 5.1f};
  std::vector<float> w{1.f, 1.f, 1.f, 1.f};
  std::vector<float> sigmaNoise{0.1f, 0.1f, 0.1f, 0.1f};
  std::vector<int> clusterIndex(x.size(), -1);

  auto queue = clue::get_queue(0u);

  // Clusterer<2>(dc, rhoc) — outlier/seeding distances optional.
  const float dc = 1.3f;
  const float rhoc = 9.0f;
  auto clusterer = clue::Clusterer<2>(dc, rhoc);

  // PointsHost<2> from per-dimension buffers (+ weights + cluster-index out buffer).
  auto points = clue::PointsHost<2>(queue, static_cast<int32_t>(x.size()), x, y, w, clusterIndex);
  points.set_density_uncertainty(std::span<float>(sigmaNoise));
  clusterer.make_clusters(points);

  const auto nclusters = points.n_clusters();
  std::span<const int32_t> seeds = clusterer.getSeeds();

  // Factory + accessors used by getClusters() in the algo.
  auto cpoints = clue::make_clustered_points<2>(
      queue, std::span<const float>(x), std::span<const float>(y), std::span<const float>(w),
      std::span<const int>(clusterIndex));
  const auto npoints = cpoints.size();
  auto idxs = cpoints.clusterIndexes();
  auto ws = cpoints.weights();
  auto c0 = cpoints.coords(0);
  auto c1 = cpoints.coords(1);

  // get_clusters -> AssociationMap (count / operator[]); centroid.
  auto clusters = clue::get_clusters(cpoints);
  const auto ncl = clusters.size();
  std::size_t hitsInFirst = (ncl > 0) ? clusters.count(0) : 0;
  if (ncl > 0) {
    auto centroid = clue::weighted_cluster_centroid<2>(cpoints, 0);
    std::cout << "centroid0 = (" << centroid[0] << ", " << centroid[1] << ")\n";
  }

  std::cout << "n_clusters=" << nclusters << " seeds=" << seeds.size() << " npoints=" << npoints
            << " assocClusters=" << ncl << " hitsInFirst=" << hitsInFirst << " ws=" << ws.size()
            << " c0=" << c0.size() << " c1=" << c1.size() << " idxs=" << idxs.size() << "\n";
  return 0;
}
