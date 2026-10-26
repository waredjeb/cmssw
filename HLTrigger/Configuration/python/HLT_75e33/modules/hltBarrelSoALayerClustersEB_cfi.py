import FWCore.ParameterSet.Config as cms

hltBarrelSoALayerClustersEB = cms.EDProducer("BarrelSoALayerClustersProducerECAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAEB"),
    barrelRecHitsLayerClustersSoA = cms.InputTag("hltBarrelCLUEsteringEB")
)
