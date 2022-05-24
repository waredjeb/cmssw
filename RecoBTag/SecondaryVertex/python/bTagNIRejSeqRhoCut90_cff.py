import FWCore.ParameterSet.Config as cms

from RecoBTag.SoftLepton.softLepton_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi import *
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import *
from RecoBTag.Combined.pfDeepCSVTagInfos_cfi import pfDeepCSVTagInfos

from RecoBTag.SecondaryVertex.candidateCombinedSecondaryVertexV2Computer_cfi import *
from RecoBTag.SecondaryVertex.pfCombinedInclusiveSecondaryVertexV2BJetTags_cfi import *
from RecoBTag.SecondaryVertex.combinedInclusiveSecondaryVertexV2BJetTags_cfi import *
from RecoBTag.Combined.pfDeepCSVJetTags_cfi import pfDeepCSVJetTags

#======================================================================
# run nuclear interaction identifications with a MODIFIED
# pfInclusiveSecondaryVertexFinderTagInfos module, so relaxing
# the cut on the SV position to rho=9.0
#======================================================================


# standard CSVIVFv2 sequence with relaxed rho cut in pfInclusiveSecondaryVertexFinderTagInfos

pfInclusiveSecondaryVertexFinderTagInfosRho90 = pfInclusiveSecondaryVertexFinderTagInfos.clone()
pfInclusiveSecondaryVertexFinderTagInfosRho90.vertexCuts.distVal2dMax = 9.0

pfCombinedInclusiveSecondaryVertexV2BJetTagsRho90 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(
        cms.InputTag("pfImpactParameterTagInfos"),
        cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosRho90")
    )
)

pfDeepCSVTagInfosRho90 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosRho90')
)

pfDeepCSVJetTagsRho90 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceStandardRho90 = cms.Sequence(
    inclusiveCandidateVertexing *
    pfImpactParameterTagInfos *
    pfInclusiveSecondaryVertexFinderTagInfosRho90 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsRho90 *
    pfDeepCSVTagInfosRho90 *
    pfDeepCSVJetTagsRho90 *
    softPFElectronsTagInfos
)




from RecoBTag.SecondaryVertex.nuclearInteractionIdentification_cff import *

# NI rejection version 0

pfImpactParameterTagInfosCleanedRho90v0 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned0")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v0 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v0"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned0")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v0 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v0"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v0"))
)

pfDeepCSVTagInfosRho90v0 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v0')
)

pfDeepCSVJetTagsRho90v0 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v0'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v0 = cms.Sequence(
    nuclearInteractionsRemoved0 *
    pfImpactParameterTagInfosCleanedRho90v0 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v0 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v0 *
    pfDeepCSVTagInfosRho90v0 *
    pfDeepCSVJetTagsRho90v0 *
    softPFElectronsTagInfos
)

# NI rejection version 1

pfImpactParameterTagInfosCleanedRho90v1 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned1")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v1 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v1"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned1")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v1 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v1"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v1"))
)

pfDeepCSVTagInfosRho90v1 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v1')
)

pfDeepCSVJetTagsRho90v1 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v1'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v1 = cms.Sequence(
    nuclearInteractionsRemoved1 *
    pfImpactParameterTagInfosCleanedRho90v1 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v1 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v1 *
    pfDeepCSVTagInfosRho90v1 *
    pfDeepCSVJetTagsRho90v1 *
    softPFElectronsTagInfos
)

# NI rejection version 2

pfImpactParameterTagInfosCleanedRho90v2 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned2")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v2 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v2"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned2")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v2 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v2"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v2"))
)

pfDeepCSVTagInfosRho90v2 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v2')
)

pfDeepCSVJetTagsRho90v2 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v2'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)

MYbtagSequenceNIremovedRho90v2 = cms.Sequence(
    nuclearInteractionsRemoved2 *
    pfImpactParameterTagInfosCleanedRho90v2 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v2 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v2 *
    pfDeepCSVTagInfosRho90v2 *
    pfDeepCSVJetTagsRho90v2 *
    softPFElectronsTagInfos
)

# NI rejection version 3

pfImpactParameterTagInfosCleanedRho90v3 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned3")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v3 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v3"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned3")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v3 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v3"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v3"))
)

pfDeepCSVTagInfosRho90v3 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v3')
)

pfDeepCSVJetTagsRho90v3 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v3'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v3 = cms.Sequence(
    nuclearInteractionsRemoved3 *
    pfImpactParameterTagInfosCleanedRho90v3 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v3 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v3 *
    pfDeepCSVTagInfosRho90v3 *
    pfDeepCSVJetTagsRho90v3 *
    softPFElectronsTagInfos
)

# NI rejection version 4

pfImpactParameterTagInfosCleanedRho90v4 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned4")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v4 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v4"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned4")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v4 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v4"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v4"))
)

pfDeepCSVTagInfosRho90v4 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v4')
)

pfDeepCSVJetTagsRho90v4 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v4'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v4 = cms.Sequence(
    nuclearInteractionsRemoved4 *
    pfImpactParameterTagInfosCleanedRho90v4 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v4 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v4 *
    pfDeepCSVTagInfosRho90v4 *
    pfDeepCSVJetTagsRho90v4 *
    softPFElectronsTagInfos
)

# NI rejection version 5

pfImpactParameterTagInfosCleanedRho90v5 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned5")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v5 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v5"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned5")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v5 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v5"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v5"))
)

pfDeepCSVTagInfosRho90v5 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v5')
)

pfDeepCSVJetTagsRho90v5 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v5'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v5 = cms.Sequence(
    nuclearInteractionsRemoved5 *
    pfImpactParameterTagInfosCleanedRho90v5 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v5 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v5 *
    pfDeepCSVTagInfosRho90v5 *
    pfDeepCSVJetTagsRho90v5 *
    softPFElectronsTagInfos
)

# NI rejection version 6

pfImpactParameterTagInfosCleanedRho90v6 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned6")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v6 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v6"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned6")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v6 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v6"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v6"))
)

pfDeepCSVTagInfosRho90v6 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v6')
)

pfDeepCSVJetTagsRho90v6 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v6'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho90v6 = cms.Sequence(
    nuclearInteractionsRemoved6 *
    pfImpactParameterTagInfosCleanedRho90v6 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v6 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v6 *
    pfDeepCSVTagInfosRho90v6 *
    pfDeepCSVJetTagsRho90v6 *
    softPFElectronsTagInfos
)


# NI rejection version 7

pfImpactParameterTagInfosCleanedRho90v7 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned7")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v7 = pfInclusiveSecondaryVertexFinderTagInfosRho90.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho90v7"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned7")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v7 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho90v7"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v7"))
)

pfDeepCSVTagInfosRho90v7 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v7')
)

pfDeepCSVJetTagsRho90v7 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho90v7'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)

MYbtagSequenceNIremovedRho90v7 = cms.Sequence(
    nuclearInteractionsRemoved7 *
    pfImpactParameterTagInfosCleanedRho90v7 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho90v7 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v7 *
    pfDeepCSVTagInfosRho90v7 *
    pfDeepCSVJetTagsRho90v7 *
    softPFElectronsTagInfos
)

# all b-tagging sequences

niRejSeqRhoCut90 = cms.Sequence(
    MYbtagSequenceStandardRho90 *
    MYbtagSequenceNIremovedRho90v0 *
    MYbtagSequenceNIremovedRho90v1 *
    MYbtagSequenceNIremovedRho90v2 *
    MYbtagSequenceNIremovedRho90v3 *
    MYbtagSequenceNIremovedRho90v4 *
    MYbtagSequenceNIremovedRho90v5 *
    MYbtagSequenceNIremovedRho90v6 *
    MYbtagSequenceNIremovedRho90v7
)
