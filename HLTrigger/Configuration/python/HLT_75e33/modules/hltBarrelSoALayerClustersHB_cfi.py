import FWCore.ParameterSet.Config as cms

hltBarrelSoALayerClustersHB = cms.EDProducer("BarrelSoALayerClustersProducerHCAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAHB"),
    barrelRecHitsLayerClustersSoA = cms.InputTag("hltBarrelCLUEsteringHB")
)
