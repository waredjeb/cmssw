import FWCore.ParameterSet.Config as cms

vertexAndTracksCandCleaned = cms.EDProducer("CandPtrProjector",
	src= cms.InputTag('particleFlow'),
	veto = cms.InputTag('nuclearInteractionCandIdentifier')
)