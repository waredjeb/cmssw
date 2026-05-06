#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUE3D.h"

#include <iterator>
#include <unordered_map>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {


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
                                                  const HGCalSoAClustersDeviceCollection& lc,
                                                  std::vector<::ticl::Trackster>& tracksters,
                                                  std::array<ticl::TICLLayerTilesDevice::View, 96> tiles) {
    using namespace cms::alpakatools;

  }

  void PatternRecognitionByCLUE3D::fillPSetDescription(::edm::ParameterSetDescription& iDesc) {
    iDesc.add<std::vector<double>>("criticalDensity", {4., 4., 4.})->setComment("Critical density in GeV");
    iDesc.add<std::vector<double>>("criticalSelfDensity", {0.15, 0.15, 0.15})
        ->setComment("Minimum ratio of self_energy/local_density to become a seed");
    iDesc.add<std::vector<int>>("densitySiblingLayers", {3, 3, 3})
        ->setComment("Layers to consider while computing local density and searching for nearest higher");
    iDesc.add<std::vector<double>>("densityEtaPhiDistanceSqr", {0.0008, 0.0008, 0.0008})
        ->setComment("Distance in eta-phi space to consider for local density");
    iDesc.add<std::vector<double>>("densityXYDistanceSqr", {3.24, 3.24, 3.24})
        ->setComment("Distance in cm on transverse plane to consider for local density");
    iDesc.add<std::vector<double>>("kernelDensityFactor", {0.2, 0.2, 0.2})
        ->setComment("Kernel factor to be applied to other LC while computing local density");
    iDesc.add<bool>("densityOnSameLayer", false);
    iDesc.add<bool>("nearestHigherOnSameLayer", false)->setComment("Allow nearestHigher to be located on same layer");
    iDesc.add<bool>("useAbsoluteProjectiveScale", true)
        ->setComment("Express all cuts in terms of r/z*z_0{,phi} projective variables");
    iDesc.add<bool>("useClusterDimensionXY", false)
        ->setComment("Use estimated cluster radius to determine compatibility while computing local density");
    iDesc.add<bool>("rescaleDensityByZ", false)->setComment("Rescale local density by extension of Z volume explored");
    iDesc.add<std::vector<double>>("criticalEtaPhiDistance", {0.025, 0.025, 0.025})
        ->setComment("Minimal distance in eta-phi space from nearestHigher to become a seed");
    iDesc.add<std::vector<double>>("criticalXYDistance", {1.8, 1.8, 1.8})
        ->setComment("Minimal distance in cm on XY plane from nearestHigher to become a seed");
    iDesc.add<std::vector<int>>("criticalZDistanceLyr", {5, 5, 5})
        ->setComment("Minimal distance in layers along Z axis from nearestHigher to become a seed");
    iDesc.add<std::vector<double>>("outlierMultiplier", {2., 2., 2.})
        ->setComment("Minimal distance in transverse space from nearestHigher to become an outlier");
    iDesc.add<std::vector<int>>("minNumLayerCluster", {2, 2, 2})->setComment("Minimum number of layer clusters");
    iDesc.add<bool>("doPidCut", false);
    iDesc.add<double>("cutHadProb", 0.5);
    iDesc.add<bool>("computeLocalTime", false);
    iDesc.add<bool>("usePCACleaning", false)->setComment("Enable PCA cleaning algorithm");
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
