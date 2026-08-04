import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.HLT_75e33.sequences.HLTBeamSpotSequence_cfi import HLTBeamSpotSequence
from HLTrigger.Configuration.HLT_75e33.sequences.HLTItLocalRecoSequence_cfi import HLTItLocalRecoSequence
from HLTrigger.Configuration.HLT_75e33.modules.hltSiPixelRecHits_cfi import hltSiPixelRecHits
from HLTrigger.Configuration.HLT_75e33.modules.hltSiPhase2RecHits_cfi import hltSiPhase2RecHits
from Validation.RecoTrack.associators_cff import hltTPClusterProducer

from SimTracker.TrackerHitAssociation.hltSimPixelTrackProducer_cff import hltSimPixelTrackProducerPhase2
from Validation.TrackingMCTruth.hltSimPixelTrackAnalyzer_cff import hltSimPixelTrackAnalyzerPhase2

simPixelTrackPhase2PreValidation = cms.Sequence(
    HLTBeamSpotSequence +       # beamspot
    HLTItLocalRecoSequence +    # local tracker reconstruction
    hltSiPixelRecHits +         # reproduce the SiPixelRecHits
    hltSiPhase2RecHits +        # produce the OT RecHits
    hltTPClusterProducer +      # run the cluster to TrackingParticle association
    hltSimPixelTrackProducerPhase2 # produce the SimPixelTracks
)

simPixelTrackPhase2Validation = cms.Sequence(
    hltSimPixelTrackAnalyzerPhase2 # SimPixelTrack validation analyzer
)
