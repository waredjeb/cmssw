import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQM_cfg import *
from Validation.HGCalValidation.hgcalValidator_cfi import hgcalValidator as _hgcalValidator
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, associatorsInstances


hgcalValidator = _hgcalValidator.clone(
    label_tst = cms.VInputTag(*[cms.InputTag(label) for label in ticlIterLabels] + [cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")]),
    allTracksterTracksterAssociatorsLabels = cms.VInputTag( *[cms.InputTag('allTrackstersToSimTrackstersAssociationsByLCs:'+associator) for associator in associatorsInstances] ),
    allTracksterTracksterByHitsAssociatorsLabels = cms.VInputTag( *[cms.InputTag('allTrackstersToSimTrackstersAssociationsByHits:'+associator) for associator in associatorsInstances] )
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(hgcalValidator,
    label_cp_fake = "mixData:MergedCaloTruth",
    label_cp_effic = "mixData:MergedCaloTruth",
    label_scl = "mixData:MergedCaloTruth",
)

from Configuration.Eras.Modifier_phase2_hgcalV10_cff import phase2_hgcalV10
phase2_hgcalV10.toModify(hgcalValidator, totallayers_to_monitor = cms.int32(50))

from Configuration.Eras.Modifier_phase2_hgcalV16_cff import phase2_hgcalV16
phase2_hgcalV16.toModify(hgcalValidator, totallayers_to_monitor = cms.int32(47))

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

lcInputMask_v5  = ["ticlTrackstersCLUE3DHigh"]
lcInputMask_v5.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])

ticl_v5.toModify(hgcalValidator,
    LayerClustersInputMask = cms.VInputTag(lcInputMask_v5),
    ticlTrackstersMerge = cms.InputTag("ticlCandidate"),
    isticlv5 = cms.untracked.bool(True),
    mergeSimToRecoAssociator = cms.InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlCandidate"),
    mergeRecoToSimAssociator = cms.InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlCandidateToticlSimTrackstersfromCPs"),
)

from Configuration.ProcessModifiers.ticl_v3_cff import ticl_v3

lcInputMask_v3 = ["ticlTrackstersEM", "ticlTrackstersHAD", "ticlTrackstersTrk", "ticlTrackstersTrkEM"]
lcInputMask_v3.extend([cms.InputTag("ticlSimTracksters", "fromCPs"), cms.InputTag("ticlSimTracksters")])

ticl_v3.toModify(hgcalValidator,
    LayerClustersInputMask = cms.VInputTag(lcInputMask_v3),
    ticlTrackstersMerge = cms.InputTag("ticlTrackstersMergeV3"),
    isticlv5 = cms.untracked.bool(True),
    mergeSimToRecoAssociator = cms.InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersMergeV3"),
    mergeRecoToSimAssociator = cms.InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersMergeV3ToticlSimTrackstersfromCPs"),
    allTracksterTracksterAssociatorsLabels = cms.VInputTag(
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersEMToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersEM',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersEMToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersEM',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersTrkEMToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersTrkEM',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersTrkEMToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersHADToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersHADToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersTrkToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersTrk',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersTrkToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersTrk',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersMergeV3ToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersMergeV3',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersMergeV3ToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersfromCPsToticlTrackstersMergeV3'
    ),
    allTracksterTracksterByHitsAssociatorsLabels = cms.VInputTag(
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersEMToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersToticlTrackstersEM',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersEMToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersfromCPsToticlTrackstersEM',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersTrkEMToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersToticlTrackstersTrkEM',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersTrkEMToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersfromCPsToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersHADToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersHADToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersfromCPsToticlTrackstersHAD',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersTrkToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersToticlTrackstersTrk',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersTrkToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersfromCPsToticlTrackstersTrk',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersMergeV3ToticlSimTracksters',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersToticlTrackstersMergeV3',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlTrackstersMergeV3ToticlSimTrackstersfromCPs',
      'allTrackstersToSimTrackstersAssociationsByHits:ticlSimTrackstersfromCPsToticlTrackstersMergeV3'
    ),
)
