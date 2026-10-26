import FWCore.ParameterSet.Config as cms

hltBarrelPFRecHitsSoAHB = cms.EDProducer("PFRecHitSoAProducerHCAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    producers = cms.VPSet(
        cms.PSet(
            src = cms.InputTag("hltBarrelRecHitsSoAHB"),
            params = cms.ESInputTag("hltESPPFRecHitHCALParams","")
        )
    ),
    topology = cms.ESInputTag("hltESPPFRecHitHCALTopology",""),
    synchronise = cms.untracked.bool(False)
)
