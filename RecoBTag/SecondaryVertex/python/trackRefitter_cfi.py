import FWCore.ParameterSet.Config as cms


vertexRefitted = cms.EDProducer("NITrackReFitter",
      secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
      primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
      beamSpot = cms.InputTag("offlineBeamSpot")
      
)
