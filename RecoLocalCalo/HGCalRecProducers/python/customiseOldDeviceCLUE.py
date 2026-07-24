import FWCore.ParameterSet.Config as cms

# cmsDriver --customise hook: on top of --procModifiers alpaka (which wires the
# device CLUEstering chain), swap the CLUEstering clustering step back to the
# original in-CMSSW device CLUE (HGCalSoARecHitsLayerClustersProducer), keeping
# the SAME module label so the rest of the device chain is untouched. Used to
# compare the three layer-clustering engines.
def customiseOldDeviceCLUE(process):
    from RecoLocalCalo.HGCalRecProducers.hgCalSoARecHitsLayerClustersProducer_cfi import hgCalSoARecHitsLayerClustersProducer
    # NOTE: the old in-CMSSW device CLUE (HGCalLayerClustersAlgoWrapper) is
    # silicon-only (hardcoded HGCalSiliconTilesConstants) and segfaults on
    # scintillator. Swap only the silicon detectors; leave HSci (BH) as the
    # device CLUEstering it was wired to.
    for det in ['EE', 'HSi', 'HFNose']:
        label = 'hgcalCLUEstering' + det
        if hasattr(process, label):
            old = getattr(process, label)
            setattr(process, label, hgCalSoARecHitsLayerClustersProducer.clone(
                hgcalRecHitsSoA = old.hgcalRecHitsSoA.value(),
                deltac = old.deltac.value(),
                kappa = old.kappa.value(),
                outlierDeltaFactor = old.outlierDeltaFactor.value(),
            ))
    return process
