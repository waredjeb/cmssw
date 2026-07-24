import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi import hgcalLayerClusters as hgcalLayerClusters_
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters as hgcalMergeLayerClusters_

from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, HGCAL_noises, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer

hgcalLayerClustersEE = hgcalLayerClusters_.clone(
    detector = 'EE',
    recHits = "HGCalRecHit:HGCEERecHits",
    plugin = dict(
        dEdXweights = HGCalRecHit.layerWeights.value(),
        #With the introduction of 7 regional factors (6 for silicon plus 1 for scintillator),
        #we extend fcPerMip (along with noises below) so that it is guaranteed that they have 6 entries.
        fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
        sciThicknessCorrection = HGCalRecHit.sciThicknessCorrection.value(),
        deltasi_index_regemfac = HGCalRecHit.deltasi_index_regemfac.value(),
        fcPerEle = fC_per_ele,
        #Extending noises as fcPerMip, see comment above.
        noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = "SiCLUE"
    )
)

hgcalLayerClustersHSi = hgcalLayerClusters_.clone(
    detector = 'FH',
    recHits = "HGCalRecHit:HGCHEFRecHits",
    plugin = dict(
        dEdXweights = HGCalRecHit.layerWeights.value(),
        #With the introduction of 7 regional factors (6 for silicon plus 1 for scintillator),
        #we extend fcPerMip (along with noises below) so that it is guaranteed that they have 6 entries.
        fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
        sciThicknessCorrection = HGCalRecHit.sciThicknessCorrection.value(),
        deltasi_index_regemfac = HGCalRecHit.deltasi_index_regemfac.value(),
        fcPerEle = fC_per_ele,
        #Extending noises as fcPerMip, see comment above.
        noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = "SiCLUE"
    )
)

hgcalLayerClustersHSci = hgcalLayerClusters_.clone(
    detector = 'BH',
    recHits = "HGCalRecHit:HGCHEBRecHits",
    plugin = dict(
        # Scintillator tiles use (eta, phi) coordinates, so the critical/seed/outlier
        # distances need the scintillator scale (the default deltac/deltas/deltao are
        # the silicon-scale values used by the EE/FH/HFNose instances).
        # deltao = 0.063 reproduces the previous effective outlier distance
        # (outlierDeltaFactor = 2.0) x (scint critical distance = 0.0315).
        deltac = [0.0315, 0.0315, 0.0315, 0.0315],
        deltas = [0.0315, 0.0315, 0.0315, 0.0315],
        deltao = [0.063, 0.063, 0.063, 0.063],
        dEdXweights = HGCalRecHit.layerWeights.value(),
        #With the introduction of 7 regional factors (6 for silicon plus 1 for scintillator),
        #we extend fcPerMip (along with noises below) so that it is guaranteed that they have 6 entries.
        fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
        sciThicknessCorrection = HGCalRecHit.sciThicknessCorrection.value(),
        deltasi_index_regemfac = HGCalRecHit.deltasi_index_regemfac.value(),
        fcPerEle = fC_per_ele,
        #Extending noises as fcPerMip, see comment above.
        noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = "SciCLUE"
    )
)

hgcalLayerClustersHFNose = hgcalLayerClusters_.clone(
    detector = 'HFNose',
    recHits = "HGCalRecHit:HGCHFNoseRecHits",
    nHitsTime = 3,
    plugin = dict(
        dEdXweights = HGCalRecHit.layerNoseWeights.value(),
        maxNumberOfThickIndices = 3,
        fcPerMip = HGCalUncalibRecHit.HGCHFNoseConfig.fCPerMIP.value(),
        thicknessCorrection = HGCalRecHit.thicknessNoseCorrection.value(),
        fcPerEle = fC_per_ele,
        noises = HGCAL_noises.values.value(),
        noiseMip = hgchebackDigitizer.digiCfg.noise.value(),
        type = "SciCLUE"
    )
)

from Configuration.Eras.Modifier_phase2_hgcalV19_cff import phase2_hgcalV19
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import fCPerMIP_mean_V19
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import nonAgedNoises_v9_v19

#The v19 geometry adds a fourth silicon sensor category (HD 200um, type 3), so
#the silicon constants have 4 entries per section (8 regional factors plus 1
#for scintillator) and the CE-H offset moves from 3 to 4. The clones above copy
#the pre-modifier defaults at import time, so the v19 values must be set here
#explicitly.
_v19SiPlugin = dict(
    deltasi_index_regemfac = 4,
    maxNumberOfThickIndices = 8,
    thicknessCorrection = [0.75, 0.76, 0.75, 0.76, 0.85, 0.85, 0.84, 0.85],
    fcPerMip = fCPerMIP_mean_V19.value() + fCPerMIP_mean_V19.value(),
    noises = nonAgedNoises_v9_v19 + nonAgedNoises_v9_v19,
)
for _clusters in (hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci):
    phase2_hgcalV19.toModify(_clusters, plugin = dict(**_v19SiPlugin))

#####################################################################
# Alpaka (device) CLUEstering chain.
#
# When the "alpaka" ProcessModifier is active, each CPU
# hgcalLayerClusters<Det> module is replaced by a device chain that ends
# in the very same legacy products (std::vector<reco::CaloCluster> at the
# default label, the "timeLayerCluster" ValueMap and, for HFNose, the
# "InitialLayerClustersMask") emitted at the SAME module label, so that
# everything downstream (hgcalMergeLayerClusters, TICL, ...) is untouched.
#
# The chain per detector is:
#   hgcalSoARecHits<Det>       (HGCRecHit  -> HGCalSoARecHits SoA)
#   hgcalCLUEstering<Det>      (CLUE clustering on the SoA)
#   hgcalSoALayerClusters<Det> (build the CaloCluster SoA)
#   _fromSoA<Det>              (SoA -> legacy reco::CaloCluster products)
#####################################################################
from RecoLocalCalo.HGCalRecProducers.hgCalSoARecHitsProducer_cfi import hgCalSoARecHitsProducer
from RecoLocalCalo.HGCalRecProducers.hgCalCLUEsteringLayerClustersProducer_cfi import hgCalCLUEsteringLayerClustersProducer
from RecoLocalCalo.HGCalRecProducers.hgCalSoALayerClustersProducer_cfi import hgCalSoALayerClustersProducer
from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersFromSoAProducer_cfi import hgCalLayerClustersFromSoAProducer

# Silicon (EE/FH/BH) share the same energy-threshold constants as the CPU
# plugins configured above (see the plugin= dicts). Scintillator (BH) also
# uses the silicon constants: its own thickness index is out of range and the
# producer then uses a zero threshold, exactly as in the CPU algorithm.
_siFcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value()
_siNoises = HGCAL_noises.values.value() + HGCAL_noises.values.value()
_siThicknessCorrection = HGCalRecHit.thicknessCorrection.value()
_siDEdXweights = HGCalRecHit.layerWeights.value()

def _makeSoARecHits(det, recHits, **kwargs):
    return hgCalSoARecHitsProducer.clone(
        detector = det,
        recHits = recHits,
        maxNumberOfThickIndices = 6,
        fcPerMip = _siFcPerMip,
        thicknessCorrection = _siThicknessCorrection,
        noises = _siNoises,
        dEdXweights = _siDEdXweights,
        **kwargs
    )

# ---- EE (silicon, electromagnetic) ----
hgcalSoARecHitsEE = _makeSoARecHits('EE', "HGCalRecHit:HGCEERecHits")
hgcalCLUEsteringEE = hgCalCLUEsteringLayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsEE', deltac = 1.3, kappa = 9., outlierDeltaFactor = 2.)
hgcalSoALayerClustersEE = hgCalSoALayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsEE', hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringEE')
_fromSoAEE = hgCalLayerClustersFromSoAProducer.clone(
    src = 'hgcalSoALayerClustersEE',
    hgcalRecHitsSoA = 'hgcalSoARecHitsEE',
    hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringEE',
    detector = 'EE')

# ---- HSi / FH (silicon, hadronic) ----
hgcalSoARecHitsHSi = _makeSoARecHits('FH', "HGCalRecHit:HGCHEFRecHits")
hgcalCLUEsteringHSi = hgCalCLUEsteringLayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSi', deltac = 1.3, kappa = 9., outlierDeltaFactor = 2.)
hgcalSoALayerClustersHSi = hgCalSoALayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSi', hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHSi')
_fromSoAFH = hgCalLayerClustersFromSoAProducer.clone(
    src = 'hgcalSoALayerClustersHSi',
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSi',
    hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHSi',
    detector = 'FH')

# ---- HSci / BH (scintillator, hadronic) ----
hgcalSoARecHitsHSci = _makeSoARecHits('BH', "HGCalRecHit:HGCHEBRecHits")
hgcalCLUEsteringHSci = hgCalCLUEsteringLayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSci', deltac = 0.0315, kappa = 9., outlierDeltaFactor = 2.)
hgcalSoALayerClustersHSci = hgCalSoALayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSci', hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHSci')
_fromSoABH = hgCalLayerClustersFromSoAProducer.clone(
    src = 'hgcalSoALayerClustersHSci',
    hgcalRecHitsSoA = 'hgcalSoARecHitsHSci',
    hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHSci',
    detector = 'BH')

# ---- HFNose (silicon, only under the phase2_hfnose modifier) ----
hgcalSoARecHitsHFNose = hgCalSoARecHitsProducer.clone(
    detector = 'HFNose',
    recHits = "HGCalRecHit:HGCHFNoseRecHits",
    maxNumberOfThickIndices = 3,
    fcPerMip = HGCalUncalibRecHit.HGCHFNoseConfig.fCPerMIP.value(),
    thicknessCorrection = HGCalRecHit.thicknessNoseCorrection.value(),
    noises = HGCAL_noises.values.value(),
    dEdXweights = HGCalRecHit.layerNoseWeights.value())
hgcalCLUEsteringHFNose = hgCalCLUEsteringLayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHFNose', deltac = 1.3, kappa = 9., outlierDeltaFactor = 2.)
hgcalSoALayerClustersHFNose = hgCalSoALayerClustersProducer.clone(
    hgcalRecHitsSoA = 'hgcalSoARecHitsHFNose', hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHFNose')
_fromSoAHFNose = hgCalLayerClustersFromSoAProducer.clone(
    src = 'hgcalSoALayerClustersHFNose',
    hgcalRecHitsSoA = 'hgcalSoARecHitsHFNose',
    hgcalRecHitsLayerClustersSoA = 'hgcalCLUEsteringHFNose',
    detector = 'HFNose',
    nHitsTime = 3)

# v19 geometry: the silicon constants gain a fourth thickness category, so the
# SoA rechit producers must be updated exactly like the CPU plugins above.
_v19SiSoAParams = dict(
    maxNumberOfThickIndices = 8,
    thicknessCorrection = [0.75, 0.76, 0.75, 0.76, 0.85, 0.85, 0.84, 0.85],
    fcPerMip = fCPerMIP_mean_V19.value() + fCPerMIP_mean_V19.value(),
    noises = nonAgedNoises_v9_v19 + nonAgedNoises_v9_v19,
)
for _soa in (hgcalSoARecHitsEE, hgcalSoARecHitsHSi, hgcalSoARecHitsHSci):
    phase2_hgcalV19.toModify(_soa, **_v19SiSoAParams)

# The device producers that must be scheduled (in a Task) under the alpaka
# modifier. Consumed collections are auto-copied device->host by the alpaka
# framework, so the CPU converter (_fromSoA*) can read the host SoAs directly.
# These Tasks are added to hgcalLocalRecoTask under alpaka in
# RecoLocalCalo/Configuration/python/hgcalLocalReco_cff.py.
hgcalLayerClustersAlpakaTask = cms.Task(
    hgcalSoARecHitsEE,   hgcalCLUEsteringEE,   hgcalSoALayerClustersEE,
    hgcalSoARecHitsHSi,  hgcalCLUEsteringHSi,  hgcalSoALayerClustersHSi,
    hgcalSoARecHitsHSci, hgcalCLUEsteringHSci, hgcalSoALayerClustersHSci,
)
hgcalLayerClustersHFNoseAlpakaTask = cms.Task(
    hgcalSoARecHitsHFNose, hgcalCLUEsteringHFNose, hgcalSoALayerClustersHFNose,
)

# Replace the CPU producers with the SoA->legacy converters at the SAME labels.
from Configuration.ProcessModifiers.alpaka_cff import alpaka
alpaka.toReplaceWith(hgcalLayerClustersEE,     _fromSoAEE)
alpaka.toReplaceWith(hgcalLayerClustersHSi,    _fromSoAFH)
alpaka.toReplaceWith(hgcalLayerClustersHSci,   _fromSoABH)
alpaka.toReplaceWith(hgcalLayerClustersHFNose, _fromSoAHFNose)

hgcalMergeLayerClusters = hgcalMergeLayerClusters_.clone(
)

layerClusters = cms.VInputTag('hgcalLayerClustersEE', 'hgcalLayerClustersHSi', 'hgcalLayerClustersHSci', 'barrelLayerClustersEB', 'barrelLayerClustersHB')
time_layerClusters = cms.VInputTag('hgcalLayerClustersEE:timeLayerCluster', 'hgcalLayerClustersHSi:timeLayerCluster', 'hgcalLayerClustersHSci:timeLayerCluster', 'barrelLayerClustersEB:timeLayerCluster', 'barrelLayerClustersHB:timeLayerCluster')
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel
ticl_barrel.toModify(hgcalMergeLayerClusters, layerClusters = layerClusters, time_layerclusters = time_layerClusters)
