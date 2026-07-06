import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer as _lcAssocByEnergyScoreProducer
from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import scAssocByEnergyScoreProducer as _scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as _layerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as _layerClusterCaloParticleAssociationProducer

from SimCalorimetry.HGCalAssociatorProducers.SimClusterToCaloParticleAssociation_cfi import SimClusterToCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import  allTrackstersToSimTrackstersAssociationsByLCs as _allTrackstersToSimTrackstersAssociationsByLCs
from SimCalorimetry.HGCalAssociatorProducers.hitToSimClusterCaloParticleAssociator_cfi import hitToSimClusterCaloParticleAssociator as _hitToSimClusterCaloParticleAssociator

from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabelsPSet as _hltTiclIterLabelsPSet

from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cff import recHitMapProducer as _recHitMapProducer

from Validation.Configuration.hltBarrelSimValid_cff import hltBarrelRecHitMapProducer as _hltBarrelRecHitMapProducer
from Validation.Configuration.hltBarrelSimValid_cff import barrel_hits
hgcal_hits = ["hltHGCalRecHit:HGCEERecHits", "hltHGCalRecHit:HGCHEFRecHits", "hltHGCalRecHit:HGCHEBRecHits"]
hltRecHitMapProducer = _hltBarrelRecHitMapProducer.clone()

hltHGCalRecHitMapProducer = _hltBarrelRecHitMapProducer.clone(
    hits = hgcal_hits,
    hgcalOnly = True,
)
from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel
(phase2_common & ~ticl_barrel).toReplaceWith(hltRecHitMapProducer, hltHGCalRecHitMapProducer)

(phase2_common & ticl_barrel).toModify(hltRecHitMapProducer,
                                       hits = [*hgcal_hits, *barrel_hits],
                                       )

from SimCalorimetry.HGCalAssociatorProducers.AllLayerClusterToTracksterAssociatorsProducer_cfi import AllLayerClusterToTracksterAssociatorsProducer as _AllLayerClusterToTracksterAssociatorsProducer

hltAllLayerClusterToTracksterAssociations = _AllLayerClusterToTracksterAssociatorsProducer.clone(
    layer_clusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabelsPSet.labels],
        cms.InputTag("hltTiclSimTracksters"),
        cms.InputTag("hltTiclSimTracksters", "fromCPs"),
    )
)

hltAllTrackstersToSimTrackstersAssociationsByLCs = _allTrackstersToSimTrackstersAssociationsByLCs.clone(
    allLCtoTSAccoc =  cms.string("hltAllLayerClusterToTracksterAssociations"),
    layerClusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabelsPSet.labels]
    ),
    simTracksterCollections = cms.VInputTag(
      cms.InputTag('hltTiclSimTracksters'),
      cms.InputTag('hltTiclSimTracksters','fromCPs')
    ),
)

# L1-seeded reco trackster labels (index into hltMergeLayerClustersL1Seeded)
_hltL1SeededRecoLabels = [l for l in _hltTiclIterLabelsPSet.labels if l.endswith("L1Seeded")]

# Dedicated L1-seeded LC->trackster maps. Both the L1-seeded reco tracksters and
# the L1-seeded SimTracksters index into hltMergeLayerClustersL1Seeded, so a single
# layer_clusters resolves everything correctly (same LC space).
hltAllLayerClusterToTracksterAssociationsL1Seeded = _AllLayerClusterToTracksterAssociatorsProducer.clone(
    layer_clusters = cms.InputTag("hltMergeLayerClustersL1Seeded"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltL1SeededRecoLabels],
        cms.InputTag("hltTiclSimTrackstersL1Seeded"),
        cms.InputTag("hltTiclSimTrackstersL1Seeded", "fromCPs"),
    )
)

# L1-seeded byLCs reco<->sim maps: L1-seeded reco vs L1-seeded SimTrackster, all in
# the L1-seeded LC space (byLCs is only meaningful within one LC space).
hltAllTrackstersToSimTrackstersAssociationsByLCsL1Seeded = _allTrackstersToSimTrackstersAssociationsByLCs.clone(
    allLCtoTSAccoc =  cms.string("hltAllLayerClusterToTracksterAssociationsL1Seeded"),
    layerClusters = cms.InputTag("hltMergeLayerClustersL1Seeded"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltL1SeededRecoLabels]
    ),
    simTracksterCollections = cms.VInputTag(
      cms.InputTag('hltTiclSimTrackstersL1Seeded'),
      cms.InputTag('hltTiclSimTrackstersL1Seeded','fromCPs')
    ),
)

from SimCalorimetry.HGCalAssociatorProducers.AllTracksterToSimTracksterAssociatorsByHitsProducer_cfi import AllTracksterToSimTracksterAssociatorsByHitsProducer as _AllTracksterToSimTracksterAssociatorsByHitsProducer

hltHitToSimClusterCaloParticleAssociator = _hitToSimClusterCaloParticleAssociator.clone(
    hitMap = 'hltHGCalRecHitMapProducer:hgcalRecHitMap',
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection'
)

from SimCalorimetry.HGCalAssociatorProducers.AllHitToTracksterAssociatorsProducer_cfi import AllHitToTracksterAssociatorsProducer as _AllHitToTracksterAssociatorsProducer

# Layer-cluster collection each trackster label indexes into: the L1-seeded
# tracksters index into hltMergeLayerClustersL1Seeded, everything else into
# hltMergeLayerClusters. This is required so the byHits associator resolves each
# trackster's vertices against its own layer clusters.
def _lcTagForLabel(label):
    return cms.InputTag("hltMergeLayerClustersL1Seeded" if label.endswith("L1Seeded") else "hltMergeLayerClusters")

# byHits hit<->trackster maps. Includes the reco tracksters, the unseeded
# SimTracksters (SC + CP) and the L1-seeded SimTracksters (SC + CP), each with
# its own layer clusters via layerClustersByCollection.
hltAllHitToTracksterAssociations =  _AllHitToTracksterAssociatorsProducer.clone(
    hitMapTag = cms.InputTag("hltHGCalRecHitMapProducer","hgcalRecHitMap"),
    hits = cms.InputTag("hltHGCalRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    layerClusters = cms.InputTag("hltMergeLayerClusters"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabelsPSet.labels],
        cms.InputTag("hltTiclSimTracksters"),
        cms.InputTag("hltTiclSimTracksters", "fromCPs"),
        cms.InputTag("hltTiclSimTrackstersL1Seeded"),
        cms.InputTag("hltTiclSimTrackstersL1Seeded", "fromCPs"),
    ),
    layerClustersByCollection = cms.VInputTag(
        *[_lcTagForLabel(label) for label in _hltTiclIterLabelsPSet.labels],
        cms.InputTag("hltMergeLayerClusters"),
        cms.InputTag("hltMergeLayerClusters"),
        cms.InputTag("hltMergeLayerClustersL1Seeded"),
        cms.InputTag("hltMergeLayerClustersL1Seeded"),
    ),
)

# byHits reco<->sim maps. The L1-seeded SimTracksters are added on the "reco"
# side so they get associated to the unseeded SimTracksters (SC + CP): the
# resulting efficiency/response is the L1-seeding-region completeness monitor.
hltAllTrackstersToSimTrackstersAssociationsByHits = _AllTracksterToSimTracksterAssociatorsByHitsProducer.clone(
    allHitToTSAccoc = cms.string("hltAllHitToTracksterAssociations"),
    hitToCaloParticleMap = cms.InputTag("hltHitToSimClusterCaloParticleAssociator","hitToCaloParticleMap"),
    hitToSimClusterMap = cms.InputTag("hltHitToSimClusterCaloParticleAssociator","hitToSimClusterMap"),
    hits = cms.InputTag("hltHGCalRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabelsPSet.labels],
        cms.InputTag("hltTiclSimTrackstersL1Seeded"),
        cms.InputTag("hltTiclSimTrackstersL1Seeded", "fromCPs"),
    ),
    simTracksterCollections = cms.VInputTag(
      'hltTiclSimTracksters',
      'hltTiclSimTracksters:fromCPs'
    ),
)

# L1-seeded byHits reco<->sim maps: L1-seeded reco vs L1-seeded SimTrackster. Reuses
# the hit maps already produced by hltAllHitToTracksterAssociations (which resolves
# the L1-seeded reco and L1-seeded SimTracksters against hltMergeLayerClustersL1Seeded)
# and the same hit->CP/SC maps (independent of the LC collection).
hltAllTrackstersToSimTrackstersAssociationsByHitsL1Seeded = _AllTracksterToSimTracksterAssociatorsByHitsProducer.clone(
    allHitToTSAccoc = cms.string("hltAllHitToTracksterAssociations"),
    hitToCaloParticleMap = cms.InputTag("hltHitToSimClusterCaloParticleAssociator","hitToCaloParticleMap"),
    hitToSimClusterMap = cms.InputTag("hltHitToSimClusterCaloParticleAssociator","hitToSimClusterMap"),
    hits = cms.InputTag("hltHGCalRecHitMapProducer", "RefProdVectorHGCRecHitCollection"),
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltL1SeededRecoLabels]
    ),
    simTracksterCollections = cms.VInputTag(
      'hltTiclSimTrackstersL1Seeded',
      'hltTiclSimTrackstersL1Seeded:fromCPs'
    ),
)
from SimCalorimetry.HGCalAssociatorProducers.hltLCToCPAssociation_cfi import hltHGCalLCToCPAssociatorByEnergyScoreProducer, hltHGCalLayerClusterCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.hltLCToSCAssociation_cfi import hltHGCalLCToSCAssociatorByEnergyScoreProducer, hltHGCalLayerClusterSimClusterAssociation

# L1-seeded LC->CaloParticle / LC->SimCluster associations: the same (hit-based)
# energy-score producers, but resolving against the L1-seeded merged layer
# clusters. These feed the L1-seeded SimTrackster (see HLTSimTracksters_cff).
hltHGCalLayerClusterCaloParticleAssociationL1Seeded = hltHGCalLayerClusterCaloParticleAssociation.clone(
    label_lc = 'hltMergeLayerClustersL1Seeded'
)
hltHGCalLayerClusterSimClusterAssociationL1Seeded = hltHGCalLayerClusterSimClusterAssociation.clone(
    label_lcl = 'hltMergeLayerClustersL1Seeded'
)

hltHgcalAssociatorsTask = cms.Task(hltHGCalRecHitMapProducer,
                                   hltHGCalLCToCPAssociatorByEnergyScoreProducer,
                                   hltHGCalLCToSCAssociatorByEnergyScoreProducer,
                                   SimClusterToCaloParticleAssociation,
                                   hltHGCalLayerClusterCaloParticleAssociation,
                                   hltHGCalLayerClusterSimClusterAssociation,
                                   hltHGCalLayerClusterCaloParticleAssociationL1Seeded,
                                   hltHGCalLayerClusterSimClusterAssociationL1Seeded,
                                   hltAllLayerClusterToTracksterAssociations,
                                   hltAllTrackstersToSimTrackstersAssociationsByLCs,
                                   hltAllHitToTracksterAssociations,
                                   hltHitToSimClusterCaloParticleAssociator,
                                   hltAllTrackstersToSimTrackstersAssociationsByHits,
                                   hltAllLayerClusterToTracksterAssociationsL1Seeded,
                                   hltAllTrackstersToSimTrackstersAssociationsByLCsL1Seeded,
                                   hltAllTrackstersToSimTrackstersAssociationsByHitsL1Seeded
                                   )

hltHgcalPrevalidation = cms.Sequence(
    hltHGCalLCToCPAssociatorByEnergyScoreProducer *
    hltHGCalLCToSCAssociatorByEnergyScoreProducer *
    hltHGCalLayerClusterCaloParticleAssociation *
    hltHGCalLayerClusterSimClusterAssociation
)
