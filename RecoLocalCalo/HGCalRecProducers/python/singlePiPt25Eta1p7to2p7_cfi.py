import FWCore.ParameterSet.Config as cms

# Minimal single-pion gun in the HGCAL acceptance (eta 1.7-2.7), for validating
# the CLUEstering-based HGCal layer-cluster reconstruction.
generator = cms.EDProducer(
    "FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(211),
        MinPt  = cms.double(25.0),
        MaxPt  = cms.double(25.0),
        MinEta = cms.double(1.7),
        MaxEta = cms.double(2.7),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(0),
    psethack  = cms.string('single pi pt 25 eta 1.7-2.7'),
    firstRun  = cms.untracked.uint32(1),
)
