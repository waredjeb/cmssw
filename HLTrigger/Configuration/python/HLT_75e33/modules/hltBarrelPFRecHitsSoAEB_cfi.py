import FWCore.ParameterSet.Config as cms

hltBarrelPFRecHitsSoAEB = cms.EDProducer("PFRecHitSoAProducerECAL@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    producers = cms.VPSet(
        cms.PSet(
            src = cms.InputTag("hltBarrelRecHitsSoAEB"),
            params = cms.ESInputTag("hltESPPFRecHitECALParams","")
        )
    ),
    topology = cms.ESInputTag("hltESPPFRecHitECALTopology",""),
    synchronise = cms.untracked.bool(False)
)
