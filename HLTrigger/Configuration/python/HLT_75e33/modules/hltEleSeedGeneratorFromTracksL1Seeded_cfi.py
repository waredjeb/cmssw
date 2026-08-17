import FWCore.ParameterSet.Config as cms

hltEleSeedGeneratorFromTracksL1Seeded = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltEleTrackSelectedByRegionL1Seeded"),
    InputVertexCollection = cms.InputTag(""),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string(''),
        propagator = cms.string('PropagatorWithMaterial')
    ),
    # originRadius/originHalfLength are only used by the usePV / InputVertexCollection
    # branches of the producer, which are disabled here; kept for reference
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    usePV = cms.bool(False),
    useEventsWithNoVertex = cms.bool(True),
    useProtoTrackKinematics = cms.bool(False),
    sortAndFilterProtoTracks = cms.bool(False),
    TTRHBuilder = cms.string('WithTrackAngle'),
    includeFourthHit = cms.bool(True),
    removeOTRechits = cms.bool(False),
    produceComplement = cms.bool(False)
)
