import FWCore.ParameterSet.Config as cms

hltBarrelPFRecHitsSoAEBL1Seeded = cms.EDProducer("PFRecHitSoAProducerECAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    producers = cms.VPSet(
        cms.PSet(
            src = cms.InputTag("hltBarrelRecHitsSoAEBL1Seeded"),
            params = cms.ESInputTag("hltESPPFRecHitECALParams","")
        )
    ),
    topology = cms.ESInputTag("hltESPPFRecHitECALTopology",""),
    synchronise = cms.untracked.bool(False)
)
