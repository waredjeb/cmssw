import FWCore.ParameterSet.Config as cms

hltTiclLayerTileProducerL1Seeded = cms.EDProducer("TICLLayerTileProducer",
    detector = cms.string('HGCAL'),
    layer_HFNose_clusters = cms.InputTag("hgcalCaloClustersFromSoAHFNose"),
    layer_clusters = cms.InputTag("hltCaloClustersFromSoAL1Seeded"),
    mightGet = cms.optional.untracked.vstring
)
