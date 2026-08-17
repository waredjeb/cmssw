import FWCore.ParameterSet.Config as cms

from ..modules.hltEgammaHoverEL1Seeded_cfi import *
from ..modules.hltEgammaElectronPixelSeedsL1Seeded_cfi import *
from ..modules.hltEgammaPixelMatchVarsL1Seeded_cfi import *
from ..modules.hltEgammaSuperClustersToPixelMatchL1Seeded_cfi import *
from ..modules.hltElePixelHitDoubletsForTripletsL1Seeded_cfi import *
from ..modules.hltElePixelHitDoubletsL1Seeded_cfi import *
from ..modules.hltElePixelHitTripletsClusterRemoverL1Seeded_cfi import *
from ..modules.hltElePixelHitTripletsL1Seeded_cfi import *
from ..modules.hltElePixelSeedsCombinedL1Seeded_cfi import *
from ..modules.hltElePixelSeedsDoubletsL1Seeded_cfi import *
from ..modules.hltElePixelSeedsTripletsL1Seeded_cfi import *
from ..modules.hltEleSeedsTrackingRegionsL1Seeded_cfi import *
from ..modules.hltPixelLayerPairsL1Seeded_cfi import *
from ..modules.hltPixelLayerTriplets_cfi import *
from ..modules.hltMeasurementTrackerEvent_cfi import *
from ..sequences.HLTDoLocalPixelSequence_cfi import *
from ..sequences.HLTDoLocalStripSequence_cfi import *

HLTElePixelMatchL1SeededSequence = cms.Sequence(HLTDoLocalPixelSequence
    +HLTDoLocalStripSequence
    +(hltMeasurementTrackerEvent
    +hltPixelLayerTriplets
    +hltEgammaHoverEL1Seeded
    +hltEgammaSuperClustersToPixelMatchL1Seeded
    +hltEleSeedsTrackingRegionsL1Seeded
    +hltElePixelHitDoubletsForTripletsL1Seeded
    +hltElePixelHitTripletsL1Seeded
    +hltElePixelSeedsTripletsL1Seeded
    +hltElePixelHitTripletsClusterRemoverL1Seeded
    +hltPixelLayerPairsL1Seeded
    +hltElePixelHitDoubletsL1Seeded
    +hltElePixelSeedsDoubletsL1Seeded
    +hltElePixelSeedsCombinedL1Seeded
    +hltEgammaElectronPixelSeedsL1Seeded
    +hltEgammaPixelMatchVarsL1Seeded))

# Seed the pixel matching from the pixel tracks instead of from the dedicated
# regional doublet/triplet seeding. The supercluster matching itself is unchanged.
from ..modules.hltEleSeedGeneratorFromTracksL1Seeded_cfi import *
from ..modules.hltEleTrackSelectedByRegionL1Seeded_cfi import *
from ..sequences.HLTItLocalRecoSequence_cfi import *
from ..sequences.HLTOtLocalRecoSequence_cfi import *
from ..sequences.HLTPhase2PixelTracksAndVerticesSequence_cfi import *

_HLTElePixelMatchL1SeededSequenceFromPixelTracks = cms.Sequence(
    HLTItLocalRecoSequence
    +HLTOtLocalRecoSequence
    +HLTPhase2PixelTracksAndVerticesSequence
    +hltEgammaHoverEL1Seeded
    +hltEgammaSuperClustersToPixelMatchL1Seeded
    +hltEleSeedsTrackingRegionsL1Seeded
    +hltEleTrackSelectedByRegionL1Seeded
    +hltEleSeedGeneratorFromTracksL1Seeded
    +hltEgammaElectronPixelSeedsL1Seeded
    +hltEgammaPixelMatchVarsL1Seeded
)

from Configuration.ProcessModifiers.hltEgammaPixelTrackSeeding_cff import hltEgammaPixelTrackSeeding
hltEgammaPixelTrackSeeding.toReplaceWith(HLTElePixelMatchL1SeededSequence,
                                         _HLTElePixelMatchL1SeededSequenceFromPixelTracks)
hltEgammaPixelTrackSeeding.toModify(hltEgammaElectronPixelSeedsL1Seeded,
                                    initialSeeds = "hltEleSeedGeneratorFromTracksL1Seeded")
