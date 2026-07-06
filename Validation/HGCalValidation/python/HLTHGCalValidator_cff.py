import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.hgcalValidator_cfi import hgcalValidator as _hgcalValidator
from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabelsPSet as _hltTiclIterLabelsPSet

hltAssociatorInstances = []
for labelts in _hltTiclIterLabelsPSet.labels:
    for labelsts in ['hltTiclSimTracksters', 'hltTiclSimTrackstersfromCPs']:
        hltAssociatorInstances.append(labelts+'To'+labelsts)
        hltAssociatorInstances.append(labelsts+'To'+labelts)

# L1-seeded SimTracksters are validated against the UNSEEDED SimTracksters to
# monitor the L1-seeding-region completeness. This comparison is cross-LC-space,
# so only a byHits association is meaningful (byLCs is impossible) -> these pairs
# go into the byHits pool only. HGCalValidator fills byLCs/byHits independently.
hltByHitsOnlyInstances = []
for labelts in ['hltTiclSimTrackstersL1Seeded', 'hltTiclSimTrackstersL1SeededfromCPs']:
    for labelsts in ['hltTiclSimTracksters', 'hltTiclSimTrackstersfromCPs']:
        hltByHitsOnlyInstances.append(labelts+'To'+labelsts)
        hltByHitsOnlyInstances.append(labelsts+'To'+labelts)

# Default is TICLv5
hltHgcalValidator = _hgcalValidator.clone(
    LayerClustersInputMask = cms.VInputTag("hltTiclTrackstersCLUE3DHigh", "hltTiclSimTracksters:fromCPs", "hltTiclSimTracksters"),
    label_tst = cms.VInputTag(*[cms.InputTag(label) for label in _hltTiclIterLabelsPSet.labels] + [cms.InputTag("hltTiclSimTracksters", "fromCPs"), cms.InputTag("hltTiclSimTracksters"), cms.InputTag("hltTiclSimTrackstersL1Seeded", "fromCPs"), cms.InputTag("hltTiclSimTrackstersL1Seeded")]),
    allTracksterTracksterAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByLCs:'+associator) for associator in hltAssociatorInstances] ),
    allTracksterTracksterByHitsAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByHits:'+associator) for associator in hltAssociatorInstances + hltByHitsOnlyInstances] ),
    associator = cms.untracked.InputTag("hltHGCalLayerClusterCaloParticleAssociation"),
    associatorSim = cms.untracked.InputTag("hltHGCalLayerClusterSimClusterAssociation"),
    dirName = cms.string('HLT/HGCAL/HGCalValidator/'),
    hits = cms.InputTag("hltHGCalRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    hitMap = cms.InputTag("hltHGCalRecHitMapProducer","hgcalRecHitMap"),
    simTrackstersMap = cms.InputTag("hltTiclSimTracksters"),
    label_layerClustersPlots = cms.string("hltHgcalMergeLayerClusters"),
    label_lcl = cms.InputTag("hltMergeLayerClusters"),
    label_simTS = cms.InputTag("hltTiclSimTracksters"),
    label_simTSFromCP = cms.InputTag("hltTiclSimTracksters","fromCPs"),
    recoTracks = cms.InputTag("hltGeneralTracks"),
    simClustersToCaloParticlesMap = cms.InputTag("SimClusterToCaloParticleAssociation","simClusterToCaloParticleMap"),
    simTiclCandidates = cms.InputTag("hltTiclSimTracksters"),
    ticlCandidates = cms.string('hltTiclCandidate'),
    ticlTrackstersMerge = cms.InputTag("hltTiclCandidate"),
    mergeRecoToSimAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs","hltTiclCandidateTohltTiclSimTrackstersfromCPs"),
    mergeSimToRecoAssociator = cms.InputTag("hltAllTrackstersToSimTrackstersAssociationsByLCs","hltTiclSimTrackstersfromCPsTohltTiclCandidate"),
)

# --- L1-seeded validator ---------------------------------------------------
# A dedicated HGCalValidator instance (the release-standard "clone retargeted to a
# different sim reference" pattern) that validates the L1-seeded reco tracksters
# against the L1-seeded SimTrackster (best-possible reconstruction WITHIN the
# L1-seeding region) instead of the unseeded SimTrackster. All associations are in
# the L1-seeded layer-cluster space. Candidate plots are off (no L1-seeded TICL
# candidate). Its own dirName -> harvested by a dedicated post-processor.
_hltL1SeededRecoLabels = [l for l in _hltTiclIterLabelsPSet.labels if l.endswith("L1Seeded")]

hltAssociatorInstancesL1Seeded = []
for labelts in _hltL1SeededRecoLabels:
    for labelsts in ['hltTiclSimTrackstersL1Seeded', 'hltTiclSimTrackstersL1SeededfromCPs']:
        hltAssociatorInstancesL1Seeded.append(labelts+'To'+labelsts)
        hltAssociatorInstancesL1Seeded.append(labelsts+'To'+labelts)

hltHgcalValidatorL1Seeded = hltHgcalValidator.clone(
    doCandidatesPlots = cms.untracked.bool(False),
    LayerClustersInputMask = cms.VInputTag("hltTiclTrackstersCLUE3DHighL1Seeded", "hltTiclSimTrackstersL1Seeded:fromCPs", "hltTiclSimTrackstersL1Seeded"),
    label_tst = cms.VInputTag(*[cms.InputTag(label) for label in _hltL1SeededRecoLabels] + [cms.InputTag("hltTiclSimTrackstersL1Seeded", "fromCPs"), cms.InputTag("hltTiclSimTrackstersL1Seeded")]),
    allTracksterTracksterAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByLCsL1Seeded:'+associator) for associator in hltAssociatorInstancesL1Seeded] ),
    allTracksterTracksterByHitsAssociatorsLabels = cms.VInputTag( *[cms.InputTag('hltAllTrackstersToSimTrackstersAssociationsByHitsL1Seeded:'+associator) for associator in hltAssociatorInstancesL1Seeded] ),
    associator = cms.untracked.InputTag("hltHGCalLayerClusterCaloParticleAssociationL1Seeded"),
    associatorSim = cms.untracked.InputTag("hltHGCalLayerClusterSimClusterAssociationL1Seeded"),
    dirName = cms.string('HLT/HGCAL/HGCalValidatorL1Seeded/'),
    simTrackstersMap = cms.InputTag("hltTiclSimTrackstersL1Seeded"),
    label_layerClustersPlots = cms.string("hltMergeLayerClustersL1Seeded"),
    label_lcl = cms.InputTag("hltMergeLayerClustersL1Seeded"),
    label_simTS = cms.InputTag("hltTiclSimTrackstersL1Seeded"),
    label_simTSFromCP = cms.InputTag("hltTiclSimTrackstersL1Seeded","fromCPs"),
)

