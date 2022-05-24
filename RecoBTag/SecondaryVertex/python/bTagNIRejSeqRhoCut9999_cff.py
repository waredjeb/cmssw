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
# pfInclusiveSecondaryVertexFinderTagInfos module, so removing
# the cut on the SV position
#======================================================================


# standard CSVIVFv2 sequence with removed rho cut in pfInclusiveSecondaryVertexFinderTagInfos

pfInclusiveSecondaryVertexFinderTagInfosRho9999 = pfInclusiveSecondaryVertexFinderTagInfos.clone()
pfInclusiveSecondaryVertexFinderTagInfosRho9999.vertexCuts.distVal2dMax = 9999.9

pfCombinedInclusiveSecondaryVertexV2BJetTagsRho9999 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(
        cms.InputTag("pfImpactParameterTagInfos"),
        cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosRho9999")
    )
)

pfDeepCSVTagInfosRho9999 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosRho9999')
)

pfDeepCSVJetTagsRho9999 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)



MYbtagSequenceStandardRho9999 = cms.Sequence(
    inclusiveCandidateVertexing *
    pfImpactParameterTagInfos *
    pfInclusiveSecondaryVertexFinderTagInfosRho9999 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsRho9999 *
    pfDeepCSVTagInfosRho9999 *
    pfDeepCSVJetTagsRho9999 *
    softPFElectronsTagInfos
)




from RecoBTag.SecondaryVertex.nuclearInteractionIdentification_cff import *

# NI rejection version 0

pfImpactParameterTagInfosCleanedRho9999v0 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned0"),
#    maximumChiSquared = 9999.9,
#    maximumLongitudinalImpactParameter = 9999.9,
#    maximumTransverseImpactParameter= 9999.9
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v0"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned0"),
)
#pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0.vertexCuts.distSig2dMin = -9999.9
#pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0.vertexCuts.massMax = -9999.9

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v0 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v0"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0"))
    
)

pfDeepCSVTagInfosRho9999v0 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0')
)

pfDeepCSVJetTagsRho9999v0 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v0'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v0 = cms.Sequence(
    nuclearInteractionsRemoved0 *
    pfImpactParameterTagInfosCleanedRho9999v0 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v0 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v0 *
    pfDeepCSVTagInfosRho9999v0 *
    pfDeepCSVJetTagsRho9999v0 *
    softPFElectronsTagInfos
)

# NI rejection version 1

pfImpactParameterTagInfosCleanedRho9999v1 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned1")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v1 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v1"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned1")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v1 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v1"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v1"))
)

pfDeepCSVTagInfosRho9999v1 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v1')
)

pfDeepCSVJetTagsRho9999v1 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v1'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v1 = cms.Sequence(
    nuclearInteractionsRemoved1 *
    pfImpactParameterTagInfosCleanedRho9999v1 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v1 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v1 *
    pfDeepCSVTagInfosRho9999v1 *
    pfDeepCSVJetTagsRho9999v1 *
    softPFElectronsTagInfos
)

# NI rejection version 2

pfImpactParameterTagInfosCleanedRho9999v2 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned2")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v2 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v2"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned2")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v2 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v2"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v2"))
)

pfDeepCSVTagInfosRho9999v2 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v2')
)

pfDeepCSVJetTagsRho9999v2 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v2'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v2 = cms.Sequence(
    nuclearInteractionsRemoved2 *
    pfImpactParameterTagInfosCleanedRho9999v2 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v2 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v2 *
    pfDeepCSVTagInfosRho9999v2 *
    pfDeepCSVJetTagsRho9999v2 *
    softPFElectronsTagInfos
)


# NI rejection version 3

pfImpactParameterTagInfosCleanedRho9999v3 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned3")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v3 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v3"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned3")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v3 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v3"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v3"))
)

pfDeepCSVTagInfosRho9999v3 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v3')
)

pfDeepCSVJetTagsRho9999v3 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v3'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v3 = cms.Sequence(
    nuclearInteractionsRemoved3 *
    pfImpactParameterTagInfosCleanedRho9999v3 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v3 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v3 *
    pfDeepCSVTagInfosRho9999v3 *
    pfDeepCSVJetTagsRho9999v3 *
    softPFElectronsTagInfos
)

# NI rejection version 4

pfImpactParameterTagInfosCleanedRho9999v4 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned4")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v4 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v4"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned4")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v4 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v4"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v4"))
)

pfDeepCSVTagInfosRho9999v4 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v4')
)

pfDeepCSVJetTagsRho9999v4 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v4'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v4 = cms.Sequence(
    nuclearInteractionsRemoved4 *
    pfImpactParameterTagInfosCleanedRho9999v4 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v4 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v4 *
    pfDeepCSVTagInfosRho9999v4 *
    pfDeepCSVJetTagsRho9999v4 *
    softPFElectronsTagInfos
)

# NI rejection version 5

pfImpactParameterTagInfosCleanedRho9999v5 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned5")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v5 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v5"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned5")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v5 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v5"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v5"))
)

pfDeepCSVTagInfosRho9999v5 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v5')
)

pfDeepCSVJetTagsRho9999v5 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v5'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v5 = cms.Sequence(
    nuclearInteractionsRemoved5 *
    pfImpactParameterTagInfosCleanedRho9999v5 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v5 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v5 *
    pfDeepCSVTagInfosRho9999v5 *
    pfDeepCSVJetTagsRho9999v5 *
    softPFElectronsTagInfos
)

# NI rejection version 6

pfImpactParameterTagInfosCleanedRho9999v6 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned6")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v6 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v6"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned6")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v6 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v6"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v6"))
)

pfDeepCSVTagInfosRho9999v6 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v6')
)

pfDeepCSVJetTagsRho9999v6 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v6'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v6 = cms.Sequence(
    nuclearInteractionsRemoved6 *
    pfImpactParameterTagInfosCleanedRho9999v6 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v6 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v6 *
    pfDeepCSVTagInfosRho9999v6 *
    pfDeepCSVJetTagsRho9999v6 *
    softPFElectronsTagInfos
)

# NI rejection version 7

pfImpactParameterTagInfosCleanedRho9999v7 = pfImpactParameterTagInfos.clone(
    candidates = cms.InputTag("vertexAndTracksCleaned7")
)

pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v7 = pfInclusiveSecondaryVertexFinderTagInfosRho9999.clone(
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v7"),
    extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned7")
)

pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v7 = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfosCleanedRho9999v7"),
                             cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v7"))
)

pfDeepCSVTagInfosRho9999v7 = pfDeepCSVTagInfos.clone(
    svTagInfos = cms.InputTag('pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v7')
)

pfDeepCSVJetTagsRho9999v7 = pfDeepCSVJetTags.clone(
    src = cms.InputTag('pfDeepCSVTagInfosRho9999v7'),
    checkSVForDefaults = cms.bool(True),
    meanPadding = cms.bool(True),
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepCSV_PhaseI.json'),
    toAdd = cms.PSet()
)


MYbtagSequenceNIremovedRho9999v7 = cms.Sequence(
    nuclearInteractionsRemoved7 *
    pfImpactParameterTagInfosCleanedRho9999v7 *
    pfInclusiveSecondaryVertexFinderTagInfosCleanedRho9999v7 *
    pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v7 *
    pfDeepCSVTagInfosRho9999v7 *
    pfDeepCSVJetTagsRho9999v7 *
    softPFElectronsTagInfos
)

# all b-tagging sequences

niRejSeqRhoCut9999 = cms.Sequence(
    MYbtagSequenceStandardRho9999 *
    MYbtagSequenceNIremovedRho9999v0 *
    MYbtagSequenceNIremovedRho9999v1 *
    MYbtagSequenceNIremovedRho9999v2 *
    MYbtagSequenceNIremovedRho9999v3 *
    MYbtagSequenceNIremovedRho9999v4 *
    MYbtagSequenceNIremovedRho9999v5 *
    MYbtagSequenceNIremovedRho9999v6 *
    MYbtagSequenceNIremovedRho9999v7
)
