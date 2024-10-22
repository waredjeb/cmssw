import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersCluesteringHigh = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    iteration_label = "CluesteringHigh"    
)

# PATTERN RECOGNITION

ticlTrackstersCluesteringHigh = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCluesteringHigh:CluesteringHigh",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CluesteringHigh",
    patternRecognitionBy = "Cluestering",
    pluginPatternRecognitionByCluestering = dict (
        criticalDensity = [0.6, 0.6, 0.6],
        criticalEtaPhiDistance = [0.025, 0.025, 0.025],
        kernelDensityFactor = [0.2, 0.2, 0.2],
        algo_verbosity = 0,
        doPidCut = True,
        cutHadProb = 999
    )

)

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
ticl_v5.toModify(ticlTrackstersCluesteringHigh.pluginPatternRecognitionByCluestering, computeLocalTime = cms.bool(True))
ticl_v5.toModify(ticlTrackstersCluesteringHigh.pluginPatternRecognitionByCluestering, usePCACleaning = cms.bool(True))

ticlCluesteringHighStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCluesteringHigh
    ,ticlTrackstersCluesteringHigh)

