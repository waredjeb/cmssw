import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalAssociatorProducers.LCToTSAssociator_cfi import layerClusterToCLUE3DTracksterAssociation, layerClusterToTracksterMergeAssociation, layerClusterToSimTracksterAssociation, layerClusterToSimTracksterFromCPsAssociation
from SimCalorimetry.HGCalAssociatorProducers.tracksterToSimTracksterAssociatorProducer_cfi import tracksterToSimTracksterAssociatorProducer
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
from Configuration.ProcessModifiers.ticl_superclustering_mustache_ticl_cff import ticl_superclustering_mustache_ticl

tracksterSimTracksterFromCPsAssociationLinking = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
    simTracksters = cms.InputTag("ticlSimTracksters", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation")
)

tracksterSimTracksterAssociationLinking = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
    simTracksters = cms.InputTag("ticlSimTracksters"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation")
)


tracksterSimTracksterFromCPsAssociationPR = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    simTracksters = cms.InputTag("ticlSimTracksters", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToCLUE3DTracksterAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation")
)

tracksterSimTracksterAssociationPR = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    simTracksters = cms.InputTag("ticlSimTracksters"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToCLUE3DTracksterAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation")
)

tracksterSimTracksterAssociationFromCPsSuperclustering = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
    simTracksters = cms.InputTag("ticlSimTracksters", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterSuperclusteringAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation")
)

tracksterSimTracksterAssociationSuperclustering = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
    simTracksters = cms.InputTag("ticlSimTracksters"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterSuperclusteringAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation")
)

tracksterSimTracksterFromCPsAssociationLinking2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation2Hits")
)

tracksterSimTracksterAssociationLinking2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterMergeAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation2Hits")
)


tracksterSimTracksterFromCPsAssociationPR2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToCLUE3DTracksterAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation2Hits")
)

tracksterSimTracksterAssociationPR2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToCLUE3DTracksterAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation2Hits")
)

tracksterSimTracksterAssociationFromCPsSuperclustering2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits", "fromCPs"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterSuperclusteringAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterFromCPsAssociation2Hits")
)

tracksterSimTracksterAssociationSuperclustering2Hits = tracksterToSimTracksterAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
    simTracksters = cms.InputTag("ticlSimTracksters2Hits"),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksterMap = cms.InputTag("layerClusterToTracksterSuperclusteringAssociation"),
    simTracksterMap = cms.InputTag("layerClusterToSimTracksterAssociation2Hits")
)

(ticl_v5 & ticl_superclustering_mustache_ticl).toModify(
    tracksterSimTracksterAssociationFromCPsSuperclustering, tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringMustache")
).toModify(
    tracksterSimTracksterAssociationSuperclustering, tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringMustache")
)


ticl_v5.toModify(tracksterSimTracksterAssociationLinking, tracksters = cms.InputTag("ticlCandidate"))
ticl_v5.toModify(tracksterSimTracksterFromCPsAssociationLinking, tracksters = cms.InputTag("ticlCandidate"))
ticl_v5.toModify(tracksterSimTracksterFromCPsAssociationLinking2Hits, tracksters = cms.InputTag("ticlCandidate"))
ticl_v5.toModify(tracksterSimTracksterAssociationLinking2Hits, tracksters = cms.InputTag("ticlCandidate"))
