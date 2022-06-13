import FWCore.ParameterSet.Config as cms

tracksterSimTracksterAssociationLinking = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(0)
)

tracksterSimTracksterAssociationPR = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(1)
)


tracksterSimTracksterAssociationPRCP = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersMerge"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(2)
)


tracksterSimTracksterAssociationLinkingbyCLUE3D = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(0)
)

tracksterSimTracksterAssociationPRbyCLUE3D = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(1)
)


tracksterSimTracksterAssociationPRCPbyCLUE3D = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
    associator = cms.InputTag('simTracksterHitLCAssociatorByEnergyScoreProducer'),
    label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    label_simTst = cms.InputTag("ticlSimTracksters"),
    label_simTstFromCP = cms.InputTag("ticlSimTracksters", "fromCPs"),
    label_lcl = cms.InputTag("hgcalLayerClusters"),
    label_scl = cms.InputTag("mix", "MergedCaloTruth"),
    label_cp = cms.InputTag("mix","MergedCaloTruth"),
    valType = cms.int32(2)
)

