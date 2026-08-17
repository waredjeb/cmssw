import FWCore.ParameterSet.Config as cms

# NOTE: these two have to be imported with 'import *', not by name. This cff is pulled
# into Configuration/StandardSequences/Validation_cff with 'import *', and Process.extend()
# labels every module it finds in that namespace before it places any sequence. Importing
# only the sequence leaves its member modules (hltOnlineBeamSpot, the local-reco modules,
# ...) unlabelled, and placing simPixelTrackPhase2PreValidation then fails with
# "An entry in sequence ... has no label" for every step3 that loads Validation_cff --
# including the ones that never request SimPixelTrackValidation.
from HLTrigger.Configuration.HLT_75e33.sequences.HLTBeamSpotSequence_cfi import *
from HLTrigger.Configuration.HLT_75e33.sequences.HLTItLocalRecoSequence_cfi import *
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
