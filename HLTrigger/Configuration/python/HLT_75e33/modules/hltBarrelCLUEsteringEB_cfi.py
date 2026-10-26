import FWCore.ParameterSet.Config as cms

hltBarrelCLUEsteringEB = cms.EDProducer("BarrelCLUEsteringLayerClustersProducerECAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAEB"),
    noise = cms.ESInputTag("hltESPPFRecHitECALParams",""),
    deltac = cms.double(1.8 * 0.0175),
    kappa = cms.double(1),
    outlierDeltaFactor = cms.double(2.)
)
