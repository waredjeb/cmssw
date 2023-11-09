import FWCore.ParameterSet.Config as cmsD

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer

# CLUSTER FILTERING/MASKING
filteredLayerClustersCLUE3DHighEM = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSizeAndLayerRange",
    min_cluster_size = 2, # inclusive
    min_layerId_ = 0, #inclusive
    max_layerId_ = 26, #inclusive
    iteration_label = "CLUE3DHighEM"
)

filteredLayerClustersCLUE3DHighHAD = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSizeAndLayerRange",
    min_cluster_size = 2, # inclusive
    min_layerId_ = 27, #inclusive
    max_layerId_ = 48, #inclusive
    iteration_label = "CLUE3DHighHAD"
)

# PATTERN RECOGNITION

ticlTrackstersCLUE3DHighEM = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCLUE3DHighEM:CLUE3DHighEM",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CLUE3DHigh",
    patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = 0.6,
        criticalEtaPhiDistance = 0.025,
        kernelDensityFactor = 0.2,
        algo_verbosity = 0
    )
)

ticlTrackstersCLUE3DHighHAD = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCLUE3DHighHAD:CLUE3DHighHAD",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CLUE3DHigh",
    patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = 0.6,
        criticalEtaPhiDistance = 0.025,
        kernelDensityFactor = 0.2,
        algo_verbosity = 0
    )
)

ticlCLUE3DHighStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCLUE3DHighEM
    ,filteredLayerClustersCLUE3DHighHAD
    ,ticlTrackstersCLUE3DHighEM
    ,ticlTrackstersCLUE3DHighHAD)

