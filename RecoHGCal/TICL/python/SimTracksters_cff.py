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
