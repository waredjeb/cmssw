import FWCore.ParameterSet.Config as cms

# HLT customisers for the HGCal device layer-clustering engine.
#
# The Phase-2 HLT (75e33) reconstructs HGCal layer clusters through the device
# SoA chain: hltHgcalSoARecHitsProducer -> hltHgcalSoARecHitsLayerClustersProducer
# (the clustering step) -> hltHgcalSoALayerClustersProducer -> converter. Each
# module has an @alpaka (device) and a SerialSync (CPU) variant.
#
# By default the clustering step is the in-CMSSW device CLUE
# (HGCalSoARecHitsLayerClustersProducer). These customisers swap ONLY the
# clustering step's module type (keeping the label, inputs and parameters), so
# the rest of the chain is untouched.

# The clustering step and its SerialSync sibling, with the target module type
# for each engine.
_CLUE_TYPES = {
    'hltHgcalSoARecHitsLayerClustersProducer': 'HGCalSoARecHitsLayerClustersProducer@alpaka',
    'hltHgcalSoARecHitsLayerClustersProducerSerialSync': 'alpaka_serial_sync::HGCalSoARecHitsLayerClustersProducer',
}
_CLUESTERING_TYPES = {
    'hltHgcalSoARecHitsLayerClustersProducer': 'HGCalCLUEsteringLayerClustersProducer@alpaka',
    'hltHgcalSoARecHitsLayerClustersProducerSerialSync': 'alpaka_serial_sync::HGCalCLUEsteringLayerClustersProducer',
}


def _retype(process, types):
    for label, newtype in types.items():
        if hasattr(process, label):
            old = getattr(process, label)
            # Rebuild the module with the new type, carrying over every parameter
            # (hgcalRecHitsSoA, deltac, kappa, outlierDeltaFactor, the alpaka PSet, ...).
            setattr(process, label, cms.EDProducer(newtype, **old.parameters_()))
    return process


def customiseHLTforCLUEstering(process):
    """Run the HGCal HLT layer clustering with the external CLUEstering device producer."""
    return _retype(process, _CLUESTERING_TYPES)


def customiseHLTforOldDeviceCLUE(process):
    """Keep the legacy in-CMSSW device CLUE for the HGCal HLT layer clustering (menu default)."""
    return _retype(process, _CLUE_TYPES)
