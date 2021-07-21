import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal
from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer as _ticlLayerTileProducer
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersLinking = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 3, # inclusive
    algo_number = 8,
    iteration_label = "Linking"
)

# CA - PATTERN RECOGNITION

ticlTrackstersLinking = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersLinking:Linking",
    seeding_regions = "ticlSeedingGlobal",
    pluginPatternRecognitionByCA = dict(
        filter_on_categories = [2, 4], # filter muons and charged hadrons
        pid_threshold = 0.0,
        skip_layers = 3,
        min_layers_per_trackster = 10,
        min_cos_theta = 0.866, # ~30 degrees
        min_cos_pointing = 0.798, # ~ 37 degrees
        max_delta_time = -1.,
        algo_verbosity = 2,
        oneTracksterPerTrackSeed = False,
        promoteEmptyRegionToTrackster = True
    ),
    itername = "Linking",
)

ticlLinkingStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersLinking
    ,ticlTrackstersLinking)

filteredLayerClustersHFNoseLinking = filteredLayerClustersLinking.clone(
    LayerClusters = 'hgcalLayerClustersHFNose',
    LayerClustersInputMask = "hgcalLayerClustersHFNose:InitialLayerClustersMask",
    min_cluster_size = 2, # inclusive
    algo_number = 9,
    iteration_label = "Linkingn"
)

ticlTrackstersHFNoseLinking = ticlTrackstersLinking.clone(
    detector = "HFNose",
    layer_clusters = "hgcalLayerClustersHFNose",
    layer_clusters_hfnose_tiles = "ticlLayerTileHFNose",
    original_mask = "hgcalLayerClustersHFNose:InitialLayerClustersMask",
    filtered_mask = "filteredLayerClustersHFNoseLinking:Linkingn",
    seeding_regions = "ticlSeedingGlobalHFNose",
    time_layerclusters = "hgcalLayerClustersHFNose:timeLayerCluster",
    itername = "Linkingn",
    pluginPatternRecognitionByCA = dict(
        filter_on_categories = [0, 1],
        min_layers_per_trackster = 5,
        pid_threshold = 0.,
        shower_start_max_layer = 5 #inclusive
    )
)

ticlHFNoseLinkingStepTask = cms.Task(ticlSeedingGlobalHFNose
    ,filteredLayerClustersHFNoseLinking
    ,ticlTrackstersHFNoseLinking)

