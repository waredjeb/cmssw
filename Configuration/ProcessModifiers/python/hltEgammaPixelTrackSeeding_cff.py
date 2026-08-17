import FWCore.ParameterSet.Config as cms

# This modifier seeds the Phase-2 HLT electron pixel matching from the pixel tracks
# (hltPhase2PixelTracks) instead of from the dedicated regional pixel doublet/triplet
# seeding, i.e. it replaces the hltElePixelSeedsCombined* chain by a region-based
# selection of the pixel tracks converted into TrajectorySeeds.
# The supercluster pixel matching itself (ElectronNHitSeedProducer) is unchanged.
hltEgammaPixelTrackSeeding = cms.Modifier()
