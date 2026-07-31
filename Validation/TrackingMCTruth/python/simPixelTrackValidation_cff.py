import FWCore.ParameterSet.Config as cms

from SimTracker.TrackerHitAssociation.simPixelTrackProducerPhase2_cfi import simPixelTrackProducerPhase2 as _simPixelTrackProducerPhase2
from Validation.TrackingMCTruth.simPixelTrackAnalyzerPhase2_cfi import simPixelTrackAnalyzerPhase2 as _simPixelTrackAnalyzerPhase2
from Validation.RecoTrack.associators_cff import hltTPClusterProducer

# SimPixelTracks are built from the HLT pixel/OT RecHits, so the modules carry the
# hlt prefix used by the HLT_75e33 menu they read from. The producer defaults
# (hltSiPixelRecHits, hltSiPhase2RecHits, hltOnlineBeamSpot) already point there.
hltSimPixelTrackProducerPhase2 = _simPixelTrackProducerPhase2.clone()

hltSimPixelTrackAnalyzerPhase2 = _simPixelTrackAnalyzerPhase2.clone(
    simPixelTrackSrc = "hltSimPixelTrackProducerPhase2"
)

# hltTPClusterProducer supplies clusterTPAssociationSrc; the phase2_tracker
# modifier repoints it at hltSiPixelClusters/hltSiPhase2Clusters.
globalPrevalidationSimPixelTrack = cms.Sequence(
    hltTPClusterProducer
    + hltSimPixelTrackProducerPhase2
)

globalValidationSimPixelTrack = cms.Sequence(
    hltSimPixelTrackAnalyzerPhase2
)
