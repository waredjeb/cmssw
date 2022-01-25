import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersEM = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSizeAndLayerRange",
    min_cluster_size = 3, # inclusive
    #max_layerId = 30,
    algo_number = 8,
    iteration_label = "EM"
)

# CA - PATTERN RECOGNITION

ticlTrackstersEM = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersEM:EM",
    seeding_regions = "ticlSeedingGlobal",
    patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalEtaPhiDistance = 0.025
    ),
    itername = "EM"
)

ticlEMStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersEM
    ,ticlTrackstersEM)

# HFNOSE CLUSTER FILTERING/MASKING

filteredLayerClustersHFNoseEM = filteredLayerClustersEM.clone(
    LayerClusters = 'hgcalLayerClustersHFNose',
    LayerClustersInputMask = 'ticlTrackstersHFNoseTrkEM',
    min_cluster_size = 3, # inclusive
    algo_number = 9,
    iteration_label = "EMn"
)

# HFNOSE CA - PATTERN RECOGNITION

ticlTrackstersHFNoseEM = ticlTrackstersEM.clone(
    detector = "HFNose",
    layer_clusters = "hgcalLayerClustersHFNose",
    layer_clusters_hfnose_tiles = "ticlLayerTileHFNose",
    original_mask = "ticlTrackstersHFNoseTrkEM",
    filtered_mask = "filteredLayerClustersHFNoseEM:EMn",
    seeding_regions = "ticlSeedingGlobalHFNose",
    time_layerclusters = "hgcalLayerClustersHFNose:timeLayerCluster",
    itername = "EMn",
    pluginPatternRecognitionByCA = dict(
       filter_on_categories = [0, 1],
       min_layers_per_trackster = 5,
       pid_threshold = 0.,
       min_cos_pointing = 0.9845, # ~10 degrees
       shower_start_max_layer = 4 ### inclusive
    )
)

ticlHFNoseEMStepTask = cms.Task(ticlSeedingGlobalHFNose
                              ,filteredLayerClustersHFNoseEM
                              ,ticlTrackstersHFNoseEM
)
