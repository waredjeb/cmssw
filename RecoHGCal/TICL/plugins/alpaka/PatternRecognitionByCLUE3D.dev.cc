// Alpaka port of PatternRecognitionbyCLUE3D
// Original author: Marco Rovere - marco.rovere@cern.ch

#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUE3D.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ticl = ::ticl;

  PatternRecognitionByCLUE3D::PatternRecognitionByCLUE3D(const edm::ParameterSet& config)
      : PatternRecognitionAlgoBase(config),
        criticalDensity_(config.getParameter<std::vector<double>>("criticalDensity")),
        criticalSelfDensity_(config.getParameter<std::vector<double>>("criticalSelfDensity")),
        densitySiblingLayers_(config.getParameter<std::vector<int>>("densitySiblingLayers")),
        densityEtaPhiDistanceSqr_(config.getParameter<std::vector<double>>("densityEtaPhiDistanceSqr")),
        densityXYDistanceSqr_(config.getParameter<std::vector<double>>("densityXYDistanceSqr")),
        kernelDensityFactor_(config.getParameter<std::vector<double>>("kernelDensityFactor")),
        densityOnSameLayer_(config.getParameter<bool>("densityOnSameLayer")),
        nearestHigherOnSameLayer_(config.getParameter<bool>("nearestHigherOnSameLayer")),
        useAbsoluteProjectiveScale_(config.getParameter<bool>("useAbsoluteProjectiveScale")),
        useClusterDimensionXY_(config.getParameter<bool>("useClusterDimensionXY")),
        rescaleDensityByZ_(config.getParameter<bool>("rescaleDensityByZ")),
        criticalEtaPhiDistance_(config.getParameter<std::vector<double>>("criticalEtaPhiDistance")),
        criticalXYDistance_(config.getParameter<std::vector<double>>("criticalXYDistance")),
        criticalZDistanceLyr_(config.getParameter<std::vector<int>>("criticalZDistanceLyr")),
        outlierMultiplier_(config.getParameter<std::vector<double>>("outlierMultiplier")),
        minNumLayerCluster_(config.getParameter<std::vector<int>>("minNumLayerCluster")),
        doPidCut_(config.getParameter<bool>("doPidCut")),
        cutHadProb_(config.getParameter<double>("cutHadProb")),
        computeLocalTime_(config.getParameter<bool>("computeLocalTime")),
        usePCACleaning_(config.getParameter<bool>("usePCACleaning")) {}

  void PatternRecognitionByCLUE3D::makeTracksters(Queue& queue,
                                                  const HGCalSoAClustersDeviceCollection& layerClusters,
                                                  const TilesConstViewArray& tiles,
                                                  std::vector<ticl::Trackster>& tracksters) {
    const int32_t nClusters = static_cast<int32_t>(layerClusters->metadata().size());

    if (nClusters == 0) {
      return;  // No clusters to process
    }

    // TODO: Implement the CLUE3D algorithm steps:
    //
    // 1. Organize clusters by layer (similar to CPU version's ClustersOnLayer)
    //    - Create device buffers for cluster properties per layer
    //    - Copy/reorganize data from layerClusters SoA
    //
    // 2. Calculate local density for each cluster
    //    - Use tiles for efficient neighbor search
    //    - Kernel: calculateLocalDensity
    //
    // 3. Calculate distance to nearest higher-density cluster
    //    - Kernel: calculateDistanceToHigher
    //
    // 4. Find seeds (high density + large distance to higher)
    //    - Kernel: findSeeds
    //
    // 5. Assign clusters to tracksters
    //    - Follow links from each cluster to its nearest-higher
    //    - Kernel: assignClusters
    //
    // 6. Copy results to host and build Trackster objects
    //    - Similar to CLUEstering implementation
    //    - Compute trackster properties (energy, barycenter, etc.)

    // For now, just a placeholder that creates empty tracksters
    tracksters.clear();

    // TODO: Remove this placeholder when implementing kernels
    // This is just skeleton code to demonstrate the structure
  }

  void PatternRecognitionByCLUE3D::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    // Critical density parameters
    iDesc.add<std::vector<double>>("criticalDensity", {0.6, 0.6, 0.6});
    iDesc.add<std::vector<double>>("criticalSelfDensity", {0.15, 0.15, 0.15});

    // Density calculation parameters
    iDesc.add<std::vector<int>>("densitySiblingLayers", {3, 3, 3});
    iDesc.add<std::vector<double>>("densityEtaPhiDistanceSqr", {0.0008, 0.0008, 0.0008});
    iDesc.add<std::vector<double>>("densityXYDistanceSqr", {3.24, 3.24, 3.24});
    iDesc.add<std::vector<double>>("kernelDensityFactor", {0.2, 0.2, 0.2});
    iDesc.add<bool>("densityOnSameLayer", false);
    iDesc.add<bool>("nearestHigherOnSameLayer", false);
    iDesc.add<bool>("useAbsoluteProjectiveScale", true);
    iDesc.add<bool>("useClusterDimensionXY", false);
    iDesc.add<bool>("rescaleDensityByZ", false);

    // Critical distance parameters
    iDesc.add<std::vector<double>>("criticalEtaPhiDistance", {0.025, 0.025, 0.025});
    iDesc.add<std::vector<double>>("criticalXYDistance", {1.8, 1.8, 1.8});
    iDesc.add<std::vector<int>>("criticalZDistanceLyr", {5, 5, 5});

    // Outlier and filtering parameters
    iDesc.add<std::vector<double>>("outlierMultiplier", {2., 2., 2.});
    iDesc.add<std::vector<int>>("minNumLayerCluster", {2, 2, 2});
    iDesc.add<bool>("doPidCut", true);
    iDesc.add<double>("cutHadProb", 0.5);

    // Additional features
    iDesc.add<bool>("computeLocalTime", false);
    iDesc.add<bool>("usePCACleaning", false);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
