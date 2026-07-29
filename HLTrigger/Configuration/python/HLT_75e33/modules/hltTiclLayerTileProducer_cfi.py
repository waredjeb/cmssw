import FWCore.ParameterSet.Config as cms

hltTiclLayerTileProducer = cms.EDProducer("TICLLayerTileProducer",
    detector = cms.string('HGCAL'),
    layer_HFNose_clusters = cms.InputTag("hgcalCaloClustersFromSoAHFNose"),
    layer_clusters = cms.InputTag("hltCaloClustersFromSoA"),
    mightGet = cms.optional.untracked.vstring
)
