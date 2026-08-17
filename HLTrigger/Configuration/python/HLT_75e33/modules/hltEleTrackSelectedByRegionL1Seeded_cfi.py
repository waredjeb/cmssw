import FWCore.ParameterSet.Config as cms

hltEleTrackSelectedByRegionL1Seeded = cms.EDProducer("TrackSelectorByRegion",
    tracks = cms.InputTag("hltPhase2PixelTracks"),
    regions = cms.InputTag("hltEleSeedsTrackingRegionsL1Seeded"),
    produceTrackCollection = cms.bool(True),
    produceMask = cms.bool(False)
)
