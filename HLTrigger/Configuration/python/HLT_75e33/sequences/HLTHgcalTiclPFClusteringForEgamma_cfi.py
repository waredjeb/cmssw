import FWCore.ParameterSet.Config as cms

from ..modules.hltFilteredLayerClustersCLUE3DHigh_cfi import *
from ..modules.hltHgcalDigis_cfi import *
from ..modules.hltMergeLayerClusters_cfi import *
from ..sequences.HLTHgcalLayerClustersSequence_cfi import *
from ..modules.hltHGCalRecHit_cfi import *
from ..modules.hltHGCalUncalibRecHit_cfi import *
from ..modules.hltParticleFlowClusterHGCal_cfi import *
from ..modules.hltParticleFlowRecHitHGC_cfi import *
from ..modules.hltParticleFlowSuperClusterHGCal_cfi import *
from ..modules.hltTiclLayerTileProducer_cfi import *
from ..modules.hltTiclSeedingGlobal_cfi import *
from ..modules.hltTiclTrackstersCLUE3DHigh_cfi import *
# TICL-barrel
from ..sequences.HLTBarrelLayerClustersSequence_cfi import *
from ..modules.hltTiclLayerTileBarrelProducer_cfi import *
from ..modules.hltTiclTrackstersCLUE3DBarrel_cfi import *

HLTHgcalTiclPFClusteringForEgamma = cms.Sequence((hltHgcalDigis+hltHGCalUncalibRecHit+hltHGCalRecHit+hltParticleFlowRecHitHGC+HLTHgcalLayerClustersSequence+hltMergeLayerClusters+hltFilteredLayerClustersCLUE3DHigh+hltTiclSeedingGlobal+hltTiclLayerTileProducer+hltTiclTrackstersCLUE3DHigh+hltParticleFlowClusterHGCal+hltParticleFlowSuperClusterHGCal))

HLTHgcalTiclPFClusteringForEgamma_barrel = cms.Sequence((
  hltHgcalDigis+
  hltHGCalUncalibRecHit+
  hltHGCalRecHit+
  hltParticleFlowRecHitHGC+
  HLTHgcalLayerClustersSequence+
  HLTBarrelLayerClustersSequence+
  hltMergeLayerClusters+
  hltFilteredLayerClustersCLUE3DHigh+
  hltTiclSeedingGlobal+
  hltTiclLayerTileProducer+
  hltTiclTrackstersCLUE3DHigh+
  hltParticleFlowClusterHGCal+
  hltParticleFlowSuperClusterHGCal+
  hltTiclLayerTileBarrelProducer+
  hltTiclTrackstersCLUE3DBarrel
))

from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel
ticl_barrel.toReplaceWith(HLTHgcalTiclPFClusteringForEgamma, HLTHgcalTiclPFClusteringForEgamma_barrel)

