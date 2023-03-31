import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer, scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as layerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as layerClusterSimClusterAssociationProducer
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import trackingParticleRecoTrackAsssociation
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import simHitTPAssocProducer
# CA - PATTERN RECOGNITION


filteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    algo_number = 8,
    min_cluster_size = 0, # inclusive
    iteration_label = "ticlSimTracksters"
)

ticlSimTracksters = _simTrackstersProducer.clone(
		filtered_mask = "filteredLayerClustersSimTracksters:ticlSimTracksters"
)


from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(ticlSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

ticlSimTrackstersTask = cms.Task(simHitTPAssocProducer, trackingParticleRecoTrackAsssociation, lcAssocByEnergyScoreProducer, layerClusterCaloParticleAssociationProducer, scAssocByEnergyScoreProducer, layerClusterSimClusterAssociationProducer, ticlSeedingGlobal, filteredLayerClustersSimTracksters, ticlSimTracksters)