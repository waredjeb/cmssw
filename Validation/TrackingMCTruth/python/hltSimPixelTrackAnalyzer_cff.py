import FWCore.ParameterSet.Config as cms

from Validation.TrackingMCTruth.simPixelTrackAnalyzerPhase2_cfi import simPixelTrackAnalyzerPhase2
from HLTrigger.Configuration.HLT_75e33.modules.hltPhase2PixelTracksSoA_cfi import hltPhase2PixelTracksSoA

hltSimPixelTrackAnalyzerPhase2 = simPixelTrackAnalyzerPhase2.clone(
    simPixelTrackSrc=cms.InputTag("hltSimPixelTrackProducerPhase2"),
    ptmin=hltPhase2PixelTracksSoA.ptmin,
    hardCurvCut=hltPhase2PixelTracksSoA.hardCurvCut,
    minHitsPerNtuplet=hltPhase2PixelTracksSoA.minHitsPerNtuplet,
    cellZ0Cut=hltPhase2PixelTracksSoA.cellZ0Cut,
    minYsizeB1=hltPhase2PixelTracksSoA.minYsizeB1,
    minYsizeB2=hltPhase2PixelTracksSoA.minYsizeB2,
    maxDYsize12=hltPhase2PixelTracksSoA.maxDYsize12,
    maxDYsize=hltPhase2PixelTracksSoA.maxDYsize,
    maxDYPred=hltPhase2PixelTracksSoA.maxDYPred,
    geometry=hltPhase2PixelTracksSoA.geometry,
)

from Configuration.ProcessModifiers.hltPhase2LegacyTrackingPatatrackQuadsChain_cff import hltPhase2LegacyTrackingPatatrackQuads
hltPhase2LegacyTrackingPatatrackQuads.toModify(hltSimPixelTrackAnalyzerPhase2, includeOTBarrel=False)
