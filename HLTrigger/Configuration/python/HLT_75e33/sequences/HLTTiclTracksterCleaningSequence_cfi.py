import FWCore.ParameterSet.Config as cms

from ..modules.hltTiclTracksterCleaning_cfi import *

HLTTiclTracksterCleaningSequence = cms.Sequence(hltTiclTracksterCleaning)