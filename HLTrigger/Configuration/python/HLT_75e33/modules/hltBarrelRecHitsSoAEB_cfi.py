import FWCore.ParameterSet.Config as cms

hltBarrelRecHitsSoAEB = cms.EDProducer("ECALRecHitSoAProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    src = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
    synchronise = cms.untracked.bool(False)
)
