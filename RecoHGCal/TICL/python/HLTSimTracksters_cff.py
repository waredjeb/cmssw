import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from Validation.RecoTrack.associators_cff import hltTrackAssociatorByHits, tpToHLTpixelTrackAssociation
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import simHitTPAssocProducer

# CA - PATTERN RECOGNITION

hltFilteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    LayerClusters = cms.InputTag("hltMergeLayerClusters"),
    LayerClustersInputMask = cms.InputTag("hltMergeLayerClusters","InitialLayerClustersMask"),
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 0, # inclusive
    iteration_label = "hltTiclSimTracksters"
)

tpToHltGeneralTrackAssociation = tpToHLTpixelTrackAssociation.clone(
    label_tr = "hltGeneralTracks"
)

hltTiclSimTracksters = _simTrackstersProducer.clone(
    layerClusterCaloParticleAssociator = cms.InputTag("hltHGCalLayerClusterCaloParticleAssociation"),
    layerClusterSimClusterAssociator = cms.InputTag("hltHGCalLayerClusterSimClusterAssociation"),
    filtered_mask = cms.InputTag("hltFilteredLayerClustersSimTracksters","hltTiclSimTracksters"),
    layer_clusters = cms.InputTag("hltMergeLayerClusters"),
    time_layerclusters = cms.InputTag("hltMergeLayerClusters","timeLayerCluster"),
    simTrackToTPMap = cms.InputTag("simHitTPAssocProducer","simTrackToTP"),
    recoTracks = cms.InputTag("hltGeneralTracks"),
    simclusters = cms.InputTag("mix","MergedCaloTruth"),
    tpToTrack = cms.InputTag("tpToHltGeneralTrackAssociation"),
    computeLocalTime = cms.bool(True)
)

from Validation.Configuration.hltHGCalSimValid_cff import *

# L1-seeded SimTracksters: the same SimTrackster reconstruction, but built from
# the L1-seeded merged layer clusters. It represents the best-possible
# reconstruction achievable within the L1-seeding region and is compared (byHits)
# to the unseeded SimTrackster to monitor how completely the region captures the
# generated particle.
hltFilteredLayerClustersSimTrackstersL1Seeded = _filteredLayerClustersProducer.clone(
    LayerClusters = cms.InputTag("hltMergeLayerClustersL1Seeded"),
    LayerClustersInputMask = cms.InputTag("hltMergeLayerClustersL1Seeded","InitialLayerClustersMask"),
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 0, # inclusive
    iteration_label = "hltTiclSimTrackstersL1Seeded"
)

hltTiclSimTrackstersL1Seeded = _simTrackstersProducer.clone(
    layerClusterCaloParticleAssociator = cms.InputTag("hltHGCalLayerClusterCaloParticleAssociationL1Seeded"),
    layerClusterSimClusterAssociator = cms.InputTag("hltHGCalLayerClusterSimClusterAssociationL1Seeded"),
    filtered_mask = cms.InputTag("hltFilteredLayerClustersSimTrackstersL1Seeded","hltTiclSimTrackstersL1Seeded"),
    layer_clusters = cms.InputTag("hltMergeLayerClustersL1Seeded"),
    time_layerclusters = cms.InputTag("hltMergeLayerClustersL1Seeded","timeLayerCluster"),
    simTrackToTPMap = cms.InputTag("simHitTPAssocProducer","simTrackToTP"),
    recoTracks = cms.InputTag("hltGeneralTracks"),
    simclusters = cms.InputTag("mix","MergedCaloTruth"),
    tpToTrack = cms.InputTag("tpToHltGeneralTrackAssociation"),
    computeLocalTime = cms.bool(True)
)

hltTiclSimTrackstersTask = cms.Task(hltTrackAssociatorByHits,
                                    tpToHltGeneralTrackAssociation,
                                    simHitTPAssocProducer,
                                    hltHgcalAssociatorsTask,
                                    hltFilteredLayerClustersSimTracksters,
                                    hltTiclSimTracksters,
                                    hltFilteredLayerClustersSimTrackstersL1Seeded,
                                    hltTiclSimTrackstersL1Seeded)

hltTiclSimTrackstersSeq = cms.Sequence(
    hltTiclSimTrackstersTask
)
