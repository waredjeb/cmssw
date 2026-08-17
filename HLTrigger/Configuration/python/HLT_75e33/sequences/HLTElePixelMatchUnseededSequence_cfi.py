import FWCore.ParameterSet.Config as cms

from ..modules.hltEgammaCandidatesUnseeded_cfi import *
from ..modules.hltEgammaElectronPixelSeedsUnseeded_cfi import *
from ..modules.hltEgammaHoverEUnseeded_cfi import *
from ..modules.hltEgammaPixelMatchVarsUnseeded_cfi import *
from ..modules.hltEgammaSuperClustersToPixelMatchUnseeded_cfi import *
from ..modules.hltElePixelHitDoubletsForTripletsUnseeded_cfi import *
from ..modules.hltElePixelHitDoubletsUnseeded_cfi import *
from ..modules.hltElePixelHitTripletsClusterRemoverUnseeded_cfi import *
from ..modules.hltElePixelHitTripletsUnseeded_cfi import *
from ..modules.hltElePixelSeedsCombinedUnseeded_cfi import *
from ..modules.hltElePixelSeedsDoubletsUnseeded_cfi import *
from ..modules.hltElePixelSeedsTripletsUnseeded_cfi import *
from ..modules.hltEleSeedsTrackingRegionsUnseeded_cfi import *
from ..modules.hltPixelLayerPairsUnseeded_cfi import *
from ..modules.hltPixelLayerTriplets_cfi import *
from ..modules.hltMeasurementTrackerEvent_cfi import *
from ..sequences.HLTDoLocalPixelSequence_cfi import *
from ..sequences.HLTDoLocalStripSequence_cfi import *

HLTElePixelMatchUnseededSequence = cms.Sequence(HLTDoLocalPixelSequence+HLTDoLocalStripSequence+(hltEgammaCandidatesUnseeded+hltEgammaHoverEUnseeded+hltMeasurementTrackerEvent+hltPixelLayerTriplets+hltEgammaSuperClustersToPixelMatchUnseeded+hltEleSeedsTrackingRegionsUnseeded+hltElePixelHitDoubletsForTripletsUnseeded+hltElePixelHitTripletsUnseeded+hltElePixelSeedsTripletsUnseeded+hltElePixelHitTripletsClusterRemoverUnseeded+hltPixelLayerPairsUnseeded+hltElePixelHitDoubletsUnseeded+hltElePixelSeedsDoubletsUnseeded+hltElePixelSeedsCombinedUnseeded+hltEgammaElectronPixelSeedsUnseeded+hltEgammaPixelMatchVarsUnseeded))

# Seed the pixel matching from the pixel tracks instead of from the dedicated
# regional doublet/triplet seeding. The supercluster matching itself is unchanged.
from ..modules.hltEleSeedGeneratorFromTracksUnseeded_cfi import *
from ..modules.hltEleTrackSelectedByRegionUnseeded_cfi import *
from ..sequences.HLTItLocalRecoSequence_cfi import *
from ..sequences.HLTOtLocalRecoSequence_cfi import *
from ..sequences.HLTPhase2PixelTracksAndVerticesSequence_cfi import *

_HLTElePixelMatchUnseededSequenceFromPixelTracks = cms.Sequence(
    HLTItLocalRecoSequence
    +HLTOtLocalRecoSequence
    +HLTPhase2PixelTracksAndVerticesSequence
    +hltEgammaCandidatesUnseeded
    +hltEgammaHoverEUnseeded
    +hltEgammaSuperClustersToPixelMatchUnseeded
    +hltEleSeedsTrackingRegionsUnseeded
    +hltEleTrackSelectedByRegionUnseeded
    +hltEleSeedGeneratorFromTracksUnseeded
    +hltEgammaElectronPixelSeedsUnseeded
    +hltEgammaPixelMatchVarsUnseeded
)

from Configuration.ProcessModifiers.hltEgammaPixelTrackSeeding_cff import hltEgammaPixelTrackSeeding
hltEgammaPixelTrackSeeding.toReplaceWith(HLTElePixelMatchUnseededSequence,
                                         _HLTElePixelMatchUnseededSequenceFromPixelTracks)
hltEgammaPixelTrackSeeding.toModify(hltEgammaElectronPixelSeedsUnseeded,
                                    initialSeeds = "hltEleSeedGeneratorFromTracksUnseeded")
