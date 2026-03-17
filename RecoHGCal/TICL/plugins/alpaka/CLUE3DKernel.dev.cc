#include <alpaka/alpaka.hpp>

#include "DataFormats/Math/interface/deltaR.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoHGCal/TICL/plugins/alpaka/CLUE3DKernel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    using namespace cms::alpakatools;

    // Device helper: calculate squared distance in eta-phi space
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static float deltaR2_etaphi(float eta1, float phi1, float eta2, float phi2) {
      float deta = eta2 - eta1;
      float dphi = phi2 - phi1;
      // Normalize phi to [-pi, pi]
      constexpr float PI = M_PI;
      if (dphi > PI)
        dphi -= 2.0f * PI;
      if (dphi < -PI)
        dphi += 2.0f * PI;
      return deta * deta + dphi * dphi;
    }

    // Device helper: calculate squared distance in XY projective coordinates
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static float distanceSqr_XY(
        float r_over_absz1, float z1, float phi1, float r_over_absz2, float z2, float phi2) {
      float r1 = r_over_absz1 * z1;
      float r2 = r_over_absz2 * z1;  // project to z1 plane
      float dphi = phi2 - phi1;
      constexpr float PI = M_PI;
      if (dphi > PI)
        dphi -= 2.0f * PI;
      if (dphi < -PI)
        dphi += 2.0f * PI;
      return (r1 - r2) * (r1 - r2) + r2 * r2 * dphi * dphi;
    }

    // Device helper: check if cluster j is reachable from cluster i for density calculation
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static bool isReachable(float r_over_absz_i,
                                                            float z_i,
                                                            float phi_i,
                                                            float r_over_absz_j,
                                                            float z_j,
                                                            float phi_j,
                                                            float distanceSqr_threshold) {
      float r_i = r_over_absz_i * z_i;
      float r_j = r_over_absz_j * z_i;  // project to z_i plane
      float dphi = phi_j - phi_i;
      constexpr float PI = M_PI;
      if (dphi > PI)
        dphi -= 2.0f * PI;
      if (dphi < -PI)
        dphi += 2.0f * PI;
      return (r_i - r_j) * (r_i - r_j) + r_j * r_j * dphi * dphi < distanceSqr_threshold;
    }

    // Device helper: get tile index from layer, eta bin, phi bin
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static int getTileIndex(int layer, int etaBin, int phiBin, int nEtaBins, int nPhiBins) {
      return layer * nEtaBins * nPhiBins + etaBin * nPhiBins + phiBin;
    }

    // Device helper: calculate eta bin from eta value
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static int getEtaBin(float eta, float etaMin, float etaMax, int nEtaBins) {
      int bin = static_cast<int>((eta - etaMin) / (etaMax - etaMin) * nEtaBins);
      return (bin < 0) ? 0 : (bin >= nEtaBins ? nEtaBins - 1 : bin);
    }

    // Device helper: calculate phi bin from phi value
    ALPAKA_FN_ACC ALPAKA_FN_INLINE static int getPhiBin(float phi, float phiMin, float phiMax, int nPhiBins) {
      int bin = static_cast<int>((phi - phiMin) / (phiMax - phiMin) * nPhiBins);
      return (bin < 0) ? 0 : (bin >= nPhiBins ? nPhiBins - 1 : bin);
    }

    //
    // Kernel 1: Calculate Local Density
    //
    struct CalculateLocalDensityKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                    const ::ticl::CLUE3DParamsSoA::ConstView params,
                                    const HGCalTilesSoA::ConstView tiles,
                                    const float* layersPosZ,
                                    CLUE3DStateSoA::View state,
                                    int nClusters) const {
        for (int i : uniform_elements(acc, nClusters)) {
          const int layer_i = state.layer(i);
          const int algoId = state.algoId(i);
          const float eta_i = state.eta(i);
          const float phi_i = state.phi(i);
          const float r_over_absz_i = state.r_over_absz(i);
          const float z_i = state.z(i);
          const int isSilicon = state.isSilicon(i);
          const int layerClusterOriginalIdx_i = state.layerClusterOriginalIdx(i);

          // Partition sides: don't mix positive and negative z
          const int lastLayerPerSide = params.lastLayerPerSide();
          int minLayer = 0;
          int maxLayer = 2 * lastLayerPerSide - 1;
          if (layer_i < lastLayerPerSide) {
            minLayer = (layer_i - params.densitySiblingLayers(algoId) < 0)
                           ? 0
                           : layer_i - params.densitySiblingLayers(algoId);
            maxLayer = (layer_i + params.densitySiblingLayers(algoId) >= lastLayerPerSide)
                           ? lastLayerPerSide - 1
                           : layer_i + params.densitySiblingLayers(algoId);
          } else {
            minLayer = (layer_i - params.densitySiblingLayers(algoId) < lastLayerPerSide)
                           ? lastLayerPerSide
                           : layer_i - params.densitySiblingLayers(algoId);
            maxLayer = (layer_i + params.densitySiblingLayers(algoId) > maxLayer)
                           ? maxLayer
                           : layer_i + params.densitySiblingLayers(algoId);
          }

          float deltaLayersZ = fabsf(layersPosZ[maxLayer % lastLayerPerSide] - layersPosZ[minLayer % lastLayerPerSide]);

          float rho = 0.0f;

          // Loop over sibling layers
          for (int currentLayer = minLayer; currentLayer <= maxLayer; ++currentLayer) {
            const bool onSameLayer = (currentLayer == layer_i);

            // Calculate eta-phi bins to search
            const int etaBinMin = getEtaBin(eta_i, params.etaMin(), params.etaMax(), params.nEtaBins()) - 2;
            const int etaBinMax = getEtaBin(eta_i, params.etaMin(), params.etaMax(), params.nEtaBins()) + 2;
            const int phiBinMin = getPhiBin(phi_i, params.phiMin(), params.phiMax(), params.nPhiBins()) - 2;
            const int phiBinMax = getPhiBin(phi_i, params.phiMin(), params.phiMax(), params.nPhiBins()) + 2;

            // Loop over eta-phi bins
            for (int ieta = (etaBinMin < 0 ? 0 : etaBinMin);
                 ieta <= (etaBinMax >= params.nEtaBins() ? params.nEtaBins() - 1 : etaBinMax);
                 ++ieta) {
              for (int iphi_it = phiBinMin; iphi_it <= phiBinMax; ++iphi_it) {
                // Handle phi periodicity
                int iphi = ((iphi_it % params.nPhiBins() + params.nPhiBins()) % params.nPhiBins());

                // Get tile range
                int tileIdx = getTileIndex(currentLayer, ieta, iphi, params.nEtaBins(), params.nPhiBins());
                int tileStart = tiles.tileOffsets(tileIdx);
                int tileEnd = tiles.tileOffsets(tileIdx + 1);

                // Loop over clusters in tile
                for (int tilePos = tileStart; tilePos < tileEnd; ++tilePos) {
                  int otherClusterIdx = tiles.tileContent(tilePos);

                  // Skip masked clusters (marked with layer == -1)
                  if (state.layer(otherClusterIdx) == -1)
                    continue;

                  const int layer_j = state.layer(otherClusterIdx);
                  const float eta_j = state.eta(otherClusterIdx);
                  const float phi_j = state.phi(otherClusterIdx);
                  const float r_over_absz_j = state.r_over_absz(otherClusterIdx);
                  const float z_j = state.z(otherClusterIdx);
                  const float energy_j = state.energy(otherClusterIdx);
                  const int layerClusterOriginalIdx_j = state.layerClusterOriginalIdx(otherClusterIdx);

                  bool onSameCluster = (layerClusterOriginalIdx_i == layerClusterOriginalIdx_j);

                  // Skip if on same layer but density not computed on same layer (unless same cluster)
                  if (onSameLayer && !params.densityOnSameLayer() && !onSameCluster)
                    continue;

                  // Check distance
                  bool reachable = false;
                  if (params.useAbsoluteProjectiveScale()) {
                    if (params.useClusterDimensionXY()) {
                      reachable = isReachable(r_over_absz_i, z_i, phi_i, r_over_absz_j, z_i, phi_j,
                                              state.radius(i) * state.radius(i));
                    } else {
                      if (isSilicon) {
                        reachable = isReachable(r_over_absz_i, z_i, phi_i, r_over_absz_j, z_i, phi_j,
                                                params.densityXYDistanceSqr(algoId));
                      } else {
                        reachable = isReachable(r_over_absz_i, z_i, phi_i, r_over_absz_j, z_i, phi_j,
                                                state.radius(i) * state.radius(i));
                      }
                    }
                  } else {
                    reachable = (deltaR2_etaphi(eta_i, phi_i, eta_j, phi_j) < params.densityEtaPhiDistanceSqr(algoId));
                  }

                  if (reachable) {
                    float factor_same_layer_different_cluster =
                        (onSameLayer && !params.densityOnSameLayer()) ? 0.0f : 1.0f;
                    float energyToAdd =
                        (onSameCluster ? 1.0f : params.kernelDensityFactor(algoId) * factor_same_layer_different_cluster) *
                        energy_j;
                    rho += energyToAdd;
                  }
                }  // end loop over clusters in tile
              }  // end loop over phi bins
            }  // end loop over eta bins
          }  // end loop over layers

          // Rescale density by Z if requested
          if (params.rescaleDensityByZ() && deltaLayersZ > 0.0f) {
            rho /= deltaLayersZ;
          }

          // Store results
          state.rho(i) = rho;
          state.z_extension(i) = deltaLayersZ;
        }  // end loop over clusters
      }
    };

    //
    // Kernel 2: Calculate Distance to Higher Density
    //
    struct CalculateDistanceToHigherKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                    const ::ticl::CLUE3DParamsSoA::ConstView params,
                                    const HGCalTilesSoA::ConstView tiles,
                                    CLUE3DStateSoA::View state,
                                    int nClusters) const {
        constexpr float maxDelta = std::numeric_limits<float>::max();
        constexpr int maxDeltaLayer = std::numeric_limits<int>::max();

        for (int i : uniform_elements(acc, nClusters)) {
          const int layer_i = state.layer(i);
          const int algoId = state.algoId(i);
          const float eta_i = state.eta(i);
          const float phi_i = state.phi(i);
          const float r_over_absz_i = state.r_over_absz(i);
          const float z_i = state.z(i);
          const float rho_i = state.rho(i);
          const int layerClusterOriginalIdx_i = state.layerClusterOriginalIdx(i);

          // Partition sides
          const int lastLayerPerSide = params.lastLayerPerSide();
          int minLayer = 0;
          int maxLayer = 2 * lastLayerPerSide - 1;
          if (layer_i < lastLayerPerSide) {
            minLayer = (layer_i - params.densitySiblingLayers(algoId) < 0)
                           ? 0
                           : layer_i - params.densitySiblingLayers(algoId);
            maxLayer = (layer_i + params.densitySiblingLayers(algoId) >= lastLayerPerSide)
                           ? lastLayerPerSide - 1
                           : layer_i + params.densitySiblingLayers(algoId);
          } else {
            minLayer = (layer_i - params.densitySiblingLayers(algoId) < lastLayerPerSide)
                           ? lastLayerPerSide
                           : layer_i - params.densitySiblingLayers(algoId);
            maxLayer = (layer_i + params.densitySiblingLayers(algoId) > maxLayer)
                           ? maxLayer
                           : layer_i + params.densitySiblingLayers(algoId);
          }

          float minDeltaDist = maxDelta;
          int minDeltaLayer = maxDeltaLayer;
          int nearestHigher_layer = -1;
          int nearestHigher_idx = -1;

          // Loop over sibling layers
          for (int currentLayer = minLayer; currentLayer <= maxLayer; ++currentLayer) {
            if (!params.nearestHigherOnSameLayer() && (currentLayer == layer_i))
              continue;

            // Calculate eta-phi bins to search (smaller window than density)
            const int etaBinMin = getEtaBin(eta_i, params.etaMin(), params.etaMax(), params.nEtaBins()) - 1;
            const int etaBinMax = getEtaBin(eta_i, params.etaMin(), params.etaMax(), params.nEtaBins()) + 1;
            const int phiBinMin = getPhiBin(phi_i, params.phiMin(), params.phiMax(), params.nPhiBins()) - 1;
            const int phiBinMax = getPhiBin(phi_i, params.phiMin(), params.phiMax(), params.nPhiBins()) + 1;

            // Loop over eta-phi bins
            for (int ieta = (etaBinMin < 0 ? 0 : etaBinMin);
                 ieta <= (etaBinMax >= params.nEtaBins() ? params.nEtaBins() - 1 : etaBinMax);
                 ++ieta) {
              for (int iphi_it = phiBinMin; iphi_it <= phiBinMax; ++iphi_it) {
                int iphi = ((iphi_it % params.nPhiBins() + params.nPhiBins()) % params.nPhiBins());

                int tileIdx = getTileIndex(currentLayer, ieta, iphi, params.nEtaBins(), params.nPhiBins());
                int tileStart = tiles.tileOffsets(tileIdx);
                int tileEnd = tiles.tileOffsets(tileIdx + 1);

                for (int tilePos = tileStart; tilePos < tileEnd; ++tilePos) {
                  int otherClusterIdx = tiles.tileContent(tilePos);

                  if (state.layer(otherClusterIdx) == -1)
                    continue;

                  const int layer_j = state.layer(otherClusterIdx);
                  const float eta_j = state.eta(otherClusterIdx);
                  const float phi_j = state.phi(otherClusterIdx);
                  const float r_over_absz_j = state.r_over_absz(otherClusterIdx);
                  const float z_j = state.z(otherClusterIdx);
                  const float rho_j = state.rho(otherClusterIdx);
                  const int layerClusterOriginalIdx_j = state.layerClusterOriginalIdx(otherClusterIdx);

                  // Check if j has higher density (or equal density but higher index for determinism)
                  bool foundHigher = (rho_j > rho_i) || (rho_j == rho_i && layerClusterOriginalIdx_j > layerClusterOriginalIdx_i);

                  if (foundHigher) {
                    float dist_transverse = 0.0f;
                    if (params.useAbsoluteProjectiveScale()) {
                      dist_transverse = distanceSqr_XY(r_over_absz_i, z_i, phi_i, r_over_absz_j, z_i, phi_j);
                      dist_transverse = sqrtf(dist_transverse);  // take sqrt for comparison
                    } else {
                      dist_transverse = sqrtf(deltaR2_etaphi(eta_i, phi_i, eta_j, phi_j));
                    }

                    int dist_layers = abs(layer_j - layer_i);

                    // Update if this is closer
                    if (dist_transverse < minDeltaDist ||
                        (dist_transverse == minDeltaDist && dist_layers < minDeltaLayer)) {
                      minDeltaDist = dist_transverse;
                      minDeltaLayer = dist_layers;
                      nearestHigher_layer = layer_j;
                      nearestHigher_idx = otherClusterIdx;
                    }
                  }
                }  // end loop over clusters in tile
              }  // end loop over phi bins
            }  // end loop over eta bins
          }  // end loop over layers

          // Store results
          state.delta_dist(i) = minDeltaDist;
          state.delta_layer(i) = minDeltaLayer;
          state.nearestHigher_layer(i) = nearestHigher_layer;
          state.nearestHigher_idx(i) = nearestHigher_idx;
        }  // end loop over clusters
      }
    };

    //
    // Kernel 3: Find Seeds and Classify Outliers
    //
    struct FindAndAssignSeedsKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                    const ::ticl::CLUE3DParamsSoA::ConstView params,
                                    CLUE3DStateSoA::View state,
                                    int nClusters,
                                    int* nSeeds) const {
        for (int i : uniform_elements(acc, nClusters)) {
          const int algoId = state.algoId(i);
          const float delta_dist = state.delta_dist(i);
          const int delta_layer = state.delta_layer(i);
          const float rho = state.rho(i);
          const float energy = state.energy(i);
          const uint8_t isSilicon = state.isSilicon(i);
          const float radius = state.radius(i);

          // Determine critical distance threshold
          const float critical_transverse_distance = params.useAbsoluteProjectiveScale()
                                                         ? params.criticalXYDistance(algoId)
                                                         : params.criticalEtaPhiDistance(algoId);

          // Seed criteria
          bool isSeed = false;
          if (isSilicon) {
            isSeed = (delta_dist > critical_transverse_distance || delta_layer > params.criticalZDistanceLyr(algoId)) &&
                     (rho >= params.criticalDensity(algoId)) &&
                     (energy / rho > params.criticalSelfDensity(algoId));
          } else {
            // Scintillator: use cluster radius instead of fixed distance
            isSeed = (delta_dist > radius || delta_layer > params.criticalZDistanceLyr(algoId)) &&
                     (rho >= params.criticalDensity(algoId)) &&
                     (energy / rho > params.criticalSelfDensity(algoId));
          }

          // Outlier criteria
          bool isOutlier = (delta_dist > params.outlierMultiplier(algoId) * critical_transverse_distance) &&
                           (rho < params.criticalDensity(algoId));

          // Assign state
          if (isSeed) {
            // Atomically increment seed counter and get cluster index
            int clusterIdx = alpaka::atomicAdd(acc, nSeeds, 1, alpaka::hierarchy::Blocks{});
            state.clusterIndex(i) = clusterIdx;
            state.isSeed(i) = 1;
            state.isOutlier(i) = 0;
          } else if (isOutlier) {
            state.clusterIndex(i) = -1;
            state.isSeed(i) = 0;
            state.isOutlier(i) = 1;
          } else {
            // Follower: will be assigned later on host
            state.clusterIndex(i) = -1;
            state.isSeed(i) = 0;
            state.isOutlier(i) = 0;
          }
        }
      }
    };

    //
    // Kernel Interface Implementation
    //

    void CLUE3DKernel::calculateLocalDensity(Queue& queue,
                                              const ::ticl::CLUE3DParamsSoA::ConstView params,
                                              const HGCalTilesSoA::ConstView tiles,
                                              const float* layersPosZ,
                                              CLUE3DStateDeviceCollection& state,
                                              int nClusters) {
      auto workDiv = make_workdiv<Acc1D>(nClusters, 256);
      alpaka::exec<Acc1D>(queue, workDiv, CalculateLocalDensityKernel{}, params, tiles, layersPosZ, state.view(), nClusters);
    }

    void CLUE3DKernel::calculateDistanceToHigher(Queue& queue,
                                                  const ::ticl::CLUE3DParamsSoA::ConstView params,
                                                  const HGCalTilesSoA::ConstView tiles,
                                                  CLUE3DStateDeviceCollection& state,
                                                  int nClusters) {
      auto workDiv = make_workdiv<Acc1D>(nClusters, 256);
      alpaka::exec<Acc1D>(
          queue, workDiv, CalculateDistanceToHigherKernel{}, params, tiles, state.view(), nClusters);
    }

    void CLUE3DKernel::findAndAssignSeeds(Queue& queue,
                                           const ::ticl::CLUE3DParamsSoA::ConstView params,
                                           CLUE3DStateDeviceCollection& state,
                                           int nClusters,
                                           int* nSeeds) {
      auto workDiv = make_workdiv<Acc1D>(nClusters, 256);
      alpaka::exec<Acc1D>(queue, workDiv, FindAndAssignSeedsKernel{}, params, state.view(), nClusters, nSeeds);
    }

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
