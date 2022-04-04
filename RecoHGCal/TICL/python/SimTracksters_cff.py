import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.slimSimTrackstersProducer_cfi import slimSimTrackstersProducer as _slimSimTrackstersProducer


# CA - PATTERN RECOGNITION

filteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    algo_number = 8,
    min_cluster_size =  0, # inclusive
    iteration_label = "ticlSimTracksters"
)

ticlSimTracksters = _simTrackstersProducer.clone(
)

ticlSlimSimTracksters = _slimSimTrackstersProducer.clone(
  detector = "HGCAL",
  layer_clusters = "hgcalLayerClusters",
  time_layerclusters = "hgcalLayerClusters:timeLayerCluster",
  filtered_mask = "ticlSimTracksters",
  patternRecognitionHighBy = "CLUE3D",
    pluginPatternRecognitionHighByCLUE3D = dict (
        criticalDensity = 1.2,
        criticalEtaPhiDistance = 0.025,
        criticalXYDistance = 3, #2.5
        minNumLayerCluster = 2,
        # densitySiblingLayers = 4,
        kernelDensityFactor = 0.2, #
        algo_verbosity = 0,
        outlierMultiplier = 10e9,
        eta_phi_window = 1
    ),
  patternRecognitionLowBy = "CLUE3D",
    pluginPatternRecognitionLowByCLUE3D = dict (
        criticalDensity = 1.5,
        criticalEtaPhiDistance = 0.025,
        criticalXYDistance = 1.8, #2.5
        minNumLayerCluster = 2,
        # densitySiblingLayers = 4,
        kernelDensityFactor = 0.2, #
        algo_verbosity = 0,
        outlierMultiplier = 10e9,
        eta_phi_window = 1
    ),
  patternRecognitionMIPBy = "CA",
        pluginPatternRecognitionMIPByCA = dict (
        skip_layers = 1,
        max_missing_layers_in_trackster = 3,
        min_layers_per_trackster = 7,
        min_cos_theta = 0.97, # ~10 degrees
        min_cos_pointing = 0.01,
        out_in_dfs = False,
        algo_verbosity = 3,
        max_delta_time = -1,
        eta_window = 1,
        phi_window = 1,
        # root_doublet_max_distance_from_seed_squared = 1, # dR=0.05,
        doSiblings = False
    )
)

## LOW 1.5GeV 1.8mm (credo)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(ticlSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

premix_stage2.toModify(ticlSlimSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

ticlSimTrackstersTask = cms.Task(filteredLayerClustersSimTracksters, ticlSimTracksters, ticlSlimSimTracksters)
