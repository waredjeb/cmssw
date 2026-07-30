import FWCore.ParameterSet.Config as cms

from ..modules.hltHgcalLayerClustersEE_cfi import *
from ..modules.hltHgcalLayerClustersHSci_cfi import *
from ..modules.hltHgcalLayerClustersHSi_cfi import *
from ..modules.hltMergeLayerClusters_cfi import *
from ..modules.hltCaloClustersFromSoA_cfi import *
from ..modules.hltHGCalRecHit_cfi import *
from ..modules.hltHGCalUncalibRecHit_cfi import *
# Heterogeneous HGCAL EE layer clusters
from ..modules.hltHgcalSoARecHitsProducer_cfi import *
from ..modules.hltHgcalSoARecHitsLayerClustersProducer_cfi import *
from ..modules.hltHgcalSoALayerClustersProducer_cfi import *
from ..modules.hltHgcalLayerClustersFromSoAProducer_cfi import *
# Heterogeneous HGCAL HSci (BH) layer clusters
from ..modules.hltHgcalSoARecHitsProducerHSci_cfi import *
from ..modules.hltHgcalSoARecHitsLayerClustersProducerHSci_cfi import *
from ..modules.hltHgcalSoALayerClustersProducerHSci_cfi import *
from ..modules.hltHgcalLayerClustersFromSoAProducerHSci_cfi import *
# Heterogeneous HGCAL HSi (FH) layer clusters
from ..modules.hltHgcalSoARecHitsProducerHSi_cfi import *
from ..modules.hltHgcalSoARecHitsLayerClustersProducerHSi_cfi import *
from ..modules.hltHgcalSoALayerClustersProducerHSi_cfi import *
from ..modules.hltHgcalLayerClustersFromSoAProducerHSi_cfi import *
# Barrel layer clusters
from ..modules.hltParticleFlowRecHitECALUnseeded_cfi import *
from ..modules.hltParticleFlowRecHitHBHE_cfi import *
from ..modules.hltBarrelLayerClustersEB_cfi import *
from ..modules.hltBarrelLayerClustersHB_cfi import *
from ..sequences.HLTPfRecHitUnseededSequence_cfi import *

from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.alpakaValidationHLT_cff import alpakaValidationHLT
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

# hltCaloClustersFromSoA kept as sequence feeds edm::Ref consumers which read AoS

HLTTICLLocalRecoSequence = cms.Sequence(
        hltHGCalUncalibRecHit+
        hltHGCalRecHit+
        hltHgcalLayerClustersEE+
        hltHgcalLayerClustersHSci+
        hltHgcalLayerClustersHSi+
        hltMergeLayerClusters+
        hltCaloClustersFromSoA)

_HLTTICLLocalRecoSequence_heterogeneous = cms.Sequence(
        hltHGCalUncalibRecHit+
        hltHGCalRecHit+
        hltHgcalSoARecHitsProducer+
        hltHgcalSoARecHitsLayerClustersProducer+
        hltHgcalSoALayerClustersProducer+
        hltHgCalLayerClustersFromSoAProducer+
        hltHgcalSoARecHitsProducerHSci+
        hltHgcalSoARecHitsLayerClustersProducerHSci+
        hltHgcalSoALayerClustersProducerHSci+
        hltHgCalLayerClustersFromSoAProducerHSci+
        hltHgcalSoARecHitsProducerHSi+
        hltHgcalSoARecHitsLayerClustersProducerHSi+
        hltHgcalSoALayerClustersProducerHSi+
        hltHgCalLayerClustersFromSoAProducerHSi+
        hltMergeLayerClusters+
        hltCaloClustersFromSoA)
(alpaka & (~ticl_barrel)).toReplaceWith(HLTTICLLocalRecoSequence, _HLTTICLLocalRecoSequence_heterogeneous)

#Define a GPU+CPU instance of TICLLocalRecoSequence, to be triggered by 'alpakaValidationHLT' procModifier
_HLTTICLLocalRecoSequence_heterogeneousGPUCPU = cms.Sequence(
        #GPU part: copied from _HLTTICLLocalRecoSequence_heterogeneous
        hltHGCalUncalibRecHit+
        hltHGCalRecHit+
        hltHgcalSoARecHitsProducer+
        hltHgcalSoARecHitsLayerClustersProducer+
        hltHgcalSoALayerClustersProducer+
        hltHgCalLayerClustersFromSoAProducer+
        hltHgcalSoARecHitsProducerHSci+
        hltHgcalSoARecHitsLayerClustersProducerHSci+
        hltHgcalSoALayerClustersProducerHSci+
        hltHgCalLayerClustersFromSoAProducerHSci+
        hltHgcalSoARecHitsProducerHSi+
        hltHgcalSoARecHitsLayerClustersProducerHSi+
        hltHgcalSoALayerClustersProducerHSi+
        hltHgCalLayerClustersFromSoAProducerHSi+
        hltHgcalLayerClustersEE+
        hltHgcalLayerClustersHSci+
        hltHgcalLayerClustersHSi+
        hltMergeLayerClusters+
        hltCaloClustersFromSoA+
        #CPU part: runs dedicated 'SerialSync' modules on CPU
        hltHgcalSoARecHitsProducerSerialSync+
        hltHgcalSoARecHitsLayerClustersProducerSerialSync+
        hltHgcalSoALayerClustersProducerSerialSync+
        hltHgCalLayerClustersFromSoAProducerSerialSync+
        hltHgcalSoARecHitsProducerHSciSerialSync+
        hltHgcalSoARecHitsLayerClustersProducerHSciSerialSync+
        hltHgcalSoALayerClustersProducerHSciSerialSync+
        hltHgCalLayerClustersFromSoAProducerHSciSerialSync+
        hltHgcalSoARecHitsProducerHSiSerialSync+
        hltHgcalSoARecHitsLayerClustersProducerHSiSerialSync+
        hltHgcalSoALayerClustersProducerHSiSerialSync+
        hltHgCalLayerClustersFromSoAProducerHSiSerialSync+
        hltMergeLayerClustersSerialSync+
        hltCaloClustersFromSoASerialSync)
alpakaValidationHLT.toReplaceWith(HLTTICLLocalRecoSequence, _HLTTICLLocalRecoSequence_heterogeneousGPUCPU)

_HLTTICLLocalRecoSequence_withBarrel = cms.Sequence(
        hltHGCalUncalibRecHit+
        hltHGCalRecHit+
        hltHgcalLayerClustersEE+
        hltHgcalLayerClustersHSci+
        hltHgcalLayerClustersHSi+
        HLTPfRecHitUnseededSequence+
        hltBarrelLayerClustersEB+
        hltBarrelLayerClustersHB+
        hltMergeLayerClusters+
        hltCaloClustersFromSoA
)
(ticl_barrel & (~alpaka)).toReplaceWith(HLTTICLLocalRecoSequence, _HLTTICLLocalRecoSequence_withBarrel)

_HLTTICLLocalRecoSequence_heterogeneous_withBarrel = cms.Sequence(
        hltHGCalUncalibRecHit+
        hltHGCalRecHit+
        hltHgcalSoARecHitsProducer+
        hltHgcalSoARecHitsLayerClustersProducer+
        hltHgcalSoALayerClustersProducer+
        hltHgCalLayerClustersFromSoAProducer+
        hltHgcalSoARecHitsProducerHSci+
        hltHgcalSoARecHitsLayerClustersProducerHSci+
        hltHgcalSoALayerClustersProducerHSci+
        hltHgCalLayerClustersFromSoAProducerHSci+
        hltHgcalSoARecHitsProducerHSi+
        hltHgcalSoARecHitsLayerClustersProducerHSi+
        hltHgcalSoALayerClustersProducerHSi+
        hltHgCalLayerClustersFromSoAProducerHSi+
        HLTPfRecHitUnseededSequence+
        hltBarrelLayerClustersEB+
        hltBarrelLayerClustersHB+
        hltMergeLayerClusters+
        hltCaloClustersFromSoA
)
(ticl_barrel & alpaka).toReplaceWith(HLTTICLLocalRecoSequence, _HLTTICLLocalRecoSequence_heterogeneous_withBarrel)
