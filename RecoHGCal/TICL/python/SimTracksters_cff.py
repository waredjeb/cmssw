import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.fineSimTrackstersProducer_cfi import fineSimTrackstersProducer as _fineSimTrackstersProducer


# CA - PATTERN RECOGNITION


filteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    algo_number = 8,
    min_cluster_size = 0, # inclusive
    iteration_label = "ticlSimTracksters"
)

ticlSimTracksters = _simTrackstersProducer.clone(
)

ticlFineSimTracksters = _fineSimTrackstersProducer.clone(
  detector = "HGCAL",
  layer_clusters = "hgcalLayerClusters",
  time_layerclusters = "hgcalLayerClusters:timeLayerCluster",
  filtered_mask = "ticlSimTracksters",
  patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = 2.,
        criticalEtaPhiDistance = 0.025,
        minNumLayerCluster = 2
    )
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(ticlSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

premix_stage2.toModify(ticlFineSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

ticlSimTrackstersTask = cms.Task(filteredLayerClustersSimTracksters, ticlSimTracksters, ticlFineSimTracksters)
# ticlSimTrackstersTask = cms.Task(filteredLayerClustersSimTracksters, ticlSimTracksters)
