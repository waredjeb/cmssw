import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cff import recHitMapProducer

# patch particle flow clusters for HGC into local reco sequence
# (for now until global reco is going with some sort of clustering)
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersHFNose, hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci, hgcalMergeLayerClusters
# Device (alpaka) CLUEstering producers. Under the alpaka modifier the
# hgcalLayerClusters<Det> modules (imported above) are already replaced, in
# hgcalLayerClusters_cff, by SoA->legacy converters that consume these
# producers; they only need to be scheduled in a Task, which is done below.
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersAlpakaTask, hgcalLayerClustersHFNoseAlpakaTask
# The individual device producers are imported as top-level names as well so
# that they are visible (via the "from ... import *" chain up to
# Reconstruction_cff) to ModuleNamesFromGlobalsVisitor when they are scheduled
# in reconstructionTask under the alpaka modifier.
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import (
    hgcalSoARecHitsEE, hgcalCLUEsteringEE, hgcalSoALayerClustersEE,
    hgcalSoARecHitsHSi, hgcalCLUEsteringHSi, hgcalSoALayerClustersHSi,
    hgcalSoARecHitsHSci, hgcalCLUEsteringHSci, hgcalSoALayerClustersHSci,
    hgcalSoARecHitsHFNose, hgcalCLUEsteringHFNose, hgcalSoALayerClustersHFNose,
)

hgcalLocalRecoTask = cms.Task( HGCalUncalibRecHit,
                                       HGCalRecHit,
                                       recHitMapProducer,
                                       hgcalLayerClustersEE,
                                       hgcalLayerClustersHSi,
                                       hgcalLayerClustersHSci,
                                       hgcalMergeLayerClusters,
                                       particleFlowRecHitHGC,
                                       particleFlowClusterHGCal )

_hfnose_hgcalLocalRecoTask = hgcalLocalRecoTask.copy()
_hfnose_hgcalLocalRecoTask.add(hgcalLayerClustersHFNose)

from Configuration.Eras.Modifier_phase2_hfnose_cff import phase2_hfnose
phase2_hfnose.toReplaceWith(
    hgcalLocalRecoTask, _hfnose_hgcalLocalRecoTask )

# Under the alpaka modifier, add the device CLUEstering producers to the local
# reco task (the hgcalLayerClusters<Det> labels are already replaced by the
# SoA->legacy converters in hgcalLayerClusters_cff).
from Configuration.ProcessModifiers.alpaka_cff import alpaka

_alpaka_hgcalLocalRecoTask = hgcalLocalRecoTask.copy()
_alpaka_hgcalLocalRecoTask.add(hgcalLayerClustersAlpakaTask)
alpaka.toReplaceWith(hgcalLocalRecoTask, _alpaka_hgcalLocalRecoTask)

# alpaka + hfnose: also schedule the HFNose device producers.
_alpaka_hfnose_hgcalLocalRecoTask = _hfnose_hgcalLocalRecoTask.copy()
_alpaka_hfnose_hgcalLocalRecoTask.add(hgcalLayerClustersAlpakaTask)
_alpaka_hfnose_hgcalLocalRecoTask.add(hgcalLayerClustersHFNoseAlpakaTask)
(phase2_hfnose & alpaka).toReplaceWith(
    hgcalLocalRecoTask, _alpaka_hfnose_hgcalLocalRecoTask )

hgcalLocalRecoSequence = cms.Sequence(hgcalLocalRecoTask)
