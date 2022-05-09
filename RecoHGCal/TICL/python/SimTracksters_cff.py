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
        criticalDensity = 0.6,
        criticalEtaPhiDistance = 0.025,
        kernelDensityFactor = 0.2,
        algo_verbosity = 0,
    ),
  patternRecognitionMIPBy = "CA",
        pluginPatternRecognitionMIPByCA = dict (
        skip_layers = 2,
        max_missing_layers_in_trackster = 3,
        min_layers_per_trackster = 7,
        min_cos_theta = 0.995, # ~10 degrees
        min_cos_pointing = 0.5,
        out_in_dfs = False,
        algo_verbosity = 3,
        max_delta_time = -1,
        # eta_window = 1,
        # phi_window = 1,
        max_out_in_hops = 1,
        # doSiblings = False
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
