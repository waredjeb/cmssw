import FWCore.ParameterSet.Config as cms

hltBarrelRecHitsSoAEBL1Seeded = cms.EDProducer("ECALRecHitSoAProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    src = cms.InputTag("hltRechitInRegionsECAL","EcalRecHitsEB"),
    synchronise = cms.untracked.bool(False)
)
