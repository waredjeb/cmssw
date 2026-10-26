import FWCore.ParameterSet.Config as cms

hltBarrelRecHitsSoAHB = cms.EDProducer("HCALRecHitSoAProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    src = cms.InputTag("hltHbhereco"),
    synchronise = cms.untracked.bool(False)
)
