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

hgcalMergeLayerClusters = hgcalMergeLayerClusters_.clone(
)

# Layer-cluster timing now travels inside the cluster SoA, so the merge no
# longer takes a parallel list of time ValueMaps.
layerClusters = cms.VInputTag('hgcalLayerClustersEE', 'hgcalLayerClustersHSi', 'hgcalLayerClustersHSci', 'barrelLayerClustersEB', 'barrelLayerClustersHB')
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel
ticl_barrel.toModify(hgcalMergeLayerClusters, layerClusters = layerClusters)

# Legacy AoS view of the merged layer clusters. TICL reads the SoA directly from
# hgcalMergeLayerClusters; this exists only for the edm::Ref consumers that cannot
# read an SoA at all (associators, validation, PF, EGamma, dumpers, Fireworks).
from RecoLocalCalo.HGCalRecProducers.caloClustersFromSoA_cfi import caloClustersFromSoA as _caloClustersFromSoA

hgcalCaloClustersFromSoA = _caloClustersFromSoA.clone(
    src = 'hgcalMergeLayerClusters'
)

# HFNose is not merged: TICL reads hgcalLayerClustersHFNose (SoA) directly. This
# converter instance exists only for the HFNose sim-truth associators.
hgcalCaloClustersFromSoAHFNose = _caloClustersFromSoA.clone(
    src = 'hgcalLayerClustersHFNose'
)
