import FWCore.ParameterSet.Config as cms

from ..modules.hltEleSeedGeneratorFromTracksL1Seeded_cfi import hltEleSeedGeneratorFromTracksL1Seeded

hltEleSeedGeneratorFromTracksUnseeded = hltEleSeedGeneratorFromTracksL1Seeded.clone(
    InputCollection = "hltEleTrackSelectedByRegionUnseeded"
)
