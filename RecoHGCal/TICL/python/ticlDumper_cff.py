import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper as ticlDumper_

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
from Configuration.ProcessModifiers.ticl_superclustering_dnn_cff import ticl_superclustering_dnn
from Configuration.ProcessModifiers.ticl_superclustering_mustache_pf_cff import ticl_superclustering_mustache_pf
from Configuration.ProcessModifiers.ticl_superclustering_mustache_ticl_cff import ticl_superclustering_mustache_ticl

ticlDumper = ticlDumper_.clone(
    saveSuperclustering = cms.bool(False),
    saveRecoSuperclusters = cms.bool(False),
    tracksterCollections = [
        cms.PSet(
            treeName=cms.string("trackstersclue3d"),
            inputTag=cms.InputTag("ticlTrackstersCLUE3DHigh")
        ),
        cms.PSet(
            treeName=cms.string("trackstersmerged"),
            inputTag=cms.InputTag("ticlTrackstersMerge")
        ),

        cms.PSet(
            treeName=cms.string("simtrackstersSC"),
            inputTag=cms.InputTag("ticlSimTracksters"),
            tracksterType=cms.string("SimTracksterSC")
        ),
        cms.PSet(
            treeName=cms.string("simtrackstersCP"),
            inputTag=cms.InputTag("ticlSimTracksters", "fromCPs"),
            tracksterType=cms.string("SimTracksterCP")
        ),
        cms.PSet(
            treeName=cms.string("simtracksters2HitsSC"),
            inputTag=cms.InputTag("ticlSimTracksters2Hits"),
            tracksterType=cms.string("SimTrackster2HitsSC")
        ),
        cms.PSet(
            treeName=cms.string("simtracksters2HitsCP"),
            inputTag=cms.InputTag("ticlSimTracksters2Hits", "fromCPs"),
            tracksterType=cms.string("SimTrackster2HitsCP")
        )
    ],
    
    associators=[
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),
        
        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("SC2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("CP2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationPR2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
        ),
        
        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("SC2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
        ),
        cms.PSet(
            branchName=cms.string("Mergetracksters"),
            suffix=cms.string("CP2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationLinking2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersMerge"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
        )
    ]
)

ticl_v5.toModify(ticlDumper,
    tracksterCollections=[
        cms.PSet(
            treeName=cms.string("trackstersCLUE3DHigh"),
            inputTag=cms.InputTag("ticlTrackstersCLUE3DHigh")
        ),
        # For future separate CLUE3D iterations :
        # cms.PSet(
        #     treeName=cms.string("trackstersCLUE3DEM"),
        #     inputTag=cms.InputTag("ticlTrackstersCLUE3DEM")
        # ),
        # cms.PSet(
        #     treeName=cms.string("trackstersCLUE3DHAD"),
        #     inputTag=cms.InputTag("ticlTrackstersCLUE3DHAD")
        # ),
        cms.PSet(
            treeName=cms.string("trackstersTiclCandidate"),
            inputTag=cms.InputTag("ticlCandidate")
        ),

        cms.PSet(
            treeName=cms.string("simtrackstersSC"),
            inputTag=cms.InputTag("ticlSimTracksters"),
            tracksterType=cms.string("SimTracksterSC")
        ),
        cms.PSet(
            treeName=cms.string("simtrackstersCP"),
            inputTag=cms.InputTag("ticlSimTracksters", "fromCPs"),
            tracksterType=cms.string("SimTracksterCP")
        ),
        cms.PSet(
            treeName=cms.string("simtracksters2HitsSC"),
            inputTag=cms.InputTag("ticlSimTracksters2Hits"),
            tracksterType=cms.string("SimTrackster2HitsSC")
        ),
        cms.PSet(
            treeName=cms.string("simtracksters2HitsCP"),
            inputTag=cms.InputTag("ticlSimTracksters2Hits", "fromCPs"),
            tracksterType=cms.string("SimTrackster2HitsCP")
        )
    ],
    ticlcandidates=cms.InputTag("ticlCandidate"),
    trackstersInCand=cms.InputTag("ticlCandidate"),
    recoSuperClusters_sourceTracksterCollection = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    associators=[
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationPR"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),

        cms.PSet(
            branchName=cms.string("ticlCandidate"),
            suffix=cms.string("SC"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlCandidate"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
        ),
        cms.PSet(
            branchName=cms.string("ticlCandidate"),
            suffix=cms.string("CP"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking"),
            tracksterCollection=cms.InputTag("ticlCandidate"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("SC2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationPR2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
        ),
        cms.PSet(
            branchName=cms.string("tsCLUE3D"),
            suffix=cms.string("CP2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationPR2Hits"),
            tracksterCollection=cms.InputTag("ticlTrackstersCLUE3DHigh"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
        ),

        cms.PSet(
            branchName=cms.string("ticlCandidate"),
            suffix=cms.string("SC2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterFromCPsAssociationLinking2Hits"),
            tracksterCollection=cms.InputTag("ticlCandidate"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
        ),
        cms.PSet(
            branchName=cms.string("ticlCandidate"),
            suffix=cms.string("CP2Hits"),
            associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationLinking2Hits"),
            tracksterCollection=cms.InputTag("ticlCandidate"),
            simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
        )
    ]
)

(ticl_v5 & ticl_superclustering_mustache_pf).toModify(ticlDumper, saveSuperclustering=False, recoSuperClusters_sourceTracksterCollection=cms.InputTag("ticlCandidate"))
(ticl_v5 & ticl_superclustering_mustache_ticl).toModify(ticlDumper, func=lambda x: x.tracksterCollections.append(cms.PSet(
            treeName=cms.string("trackstersSuperclusteringMustache"),
            inputTag=cms.InputTag("ticlTracksterLinksSuperclusteringMustache")
        ))
).toModify(ticlDumper, func=lambda x: x.associators.extend([
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("SC"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationSuperclustering"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringMustache"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("CP"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationFromCPsSuperclustering"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringMustache"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("SC"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationSuperclustering2Hits"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringMustache"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("CP"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationFromCPsSuperclustering2Hits"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringMustache"),
        simTracksterCollection=cms.InputTag("ticlSimTracksteris2Hits", "fromCPs")
    ),
])).toModify(ticlDumper,
    superclustering=cms.InputTag("ticlTracksterLinksSuperclusteringMustache")
)
(ticl_v5 & ticl_superclustering_dnn).toModify(ticlDumper, func=lambda x: x.tracksterCollections.append(cms.PSet(
            treeName=cms.string("trackstersSuperclusteringDNN"),
            inputTag=cms.InputTag("ticlTracksterLinksSuperclusteringDNN")
        ))
).toModify(ticlDumper, func=lambda x: x.associators.extend([
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("SC"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationSuperclustering"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("CP"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationFromCPsSuperclustering"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters", "fromCPs")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("SC"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationSuperclustering2Hits"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits")
    ),
    cms.PSet(
        branchName=cms.string("tsSuperclusters"),
        suffix=cms.string("CP"),
        associatorInputTag=cms.InputTag("tracksterSimTracksterAssociationFromCPsSuperclustering2Hits"),
        tracksterCollection=cms.InputTag("ticlTracksterLinksSuperclusteringDNN"),
        simTracksterCollection=cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
    ),
])).toModify(ticlDumper,
    superclustering=cms.InputTag("ticlTracksterLinksSuperclusteringDNN")
)
