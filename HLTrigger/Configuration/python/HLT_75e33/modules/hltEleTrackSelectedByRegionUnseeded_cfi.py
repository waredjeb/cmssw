import FWCore.ParameterSet.Config as cms

hltEleTrackSelectedByRegionUnseeded = cms.EDProducer("TrackSelectorByRegion",
    tracks = cms.InputTag("hltPhase2PixelTracks"),
    regions = cms.InputTag("hltEleSeedsTrackingRegionsUnseeded"),
    produceTrackCollection = cms.bool(True),
    produceMask = cms.bool(False)
)
