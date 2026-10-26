import FWCore.ParameterSet.Config as cms

hltBarrelSoALayerClustersEBL1Seeded = cms.EDProducer("BarrelSoALayerClustersProducerECAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAEBL1Seeded"),
    barrelRecHitsLayerClustersSoA = cms.InputTag("hltBarrelCLUEsteringEBL1Seeded")
)
