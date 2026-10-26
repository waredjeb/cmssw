import FWCore.ParameterSet.Config as cms

hltBarrelCLUEsteringHB = cms.EDProducer("BarrelCLUEsteringLayerClustersProducerHCAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAHB"),
    noise = cms.ESInputTag("hltESPPFRecHitHCALTopology",""),
    deltac = cms.double(3 * 0.087),
    kappa = cms.double(0),
    outlierDeltaFactor = cms.double(2.)
)
