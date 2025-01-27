import FWCore.ParameterSet.Config as cms

from ..modules.hltTiclTracksterLinksBySkeletons_cfi import *

HLTTiclTracksterLinksSequence = cms.Sequence(hltTiclTracksterLinksBySkeletonsUnseeded)
