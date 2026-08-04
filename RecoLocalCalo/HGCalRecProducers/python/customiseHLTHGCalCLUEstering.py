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


def customiseHLTforCLUEsteringAllDetectors(process):
    """Run the HGCal HLT layer clustering with CLUEstering on EE, FH (HSi) AND BH (HSci).

    The 75e33 menu only has a device SoA chain for EE (FH/BH use the CPU
    HGCalLayerClusterProducer). This customiser (1) retypes the EE clustering step
    to CLUEstering, then (2) builds full device chains for FH and BH by cloning the
    EE device modules, (3) repoints hltMergeLayerClusters at the FH/BH SoA->legacy
    converters, and (4) schedules the new producers. CLUEstering handles the
    scintillator (BH) natively (eta/phi coordinates); the old in-CMSSW device CLUE
    cannot (it is silicon-tiled).
    """
    # EE: retype the existing device clustering step to CLUEstering.
    process = customiseHLTforCLUEstering(process)

    # Scintillator (BH) noise, in MIP units. CLUEstering's effective seeding
    # density is min_density * sigmaNoise, so the BH rechit SoA must carry a
    # non-zero scintillator sigmaNoise. Take the exact menu values from the CPU
    # scintillator layer-cluster module when present, else fall back to the
    # producer defaults.
    # These parameters do not exist yet on the (menu-dumped) SoA rechit module,
    # so they must be passed to clone() as cms types, not bare Python values.
    # In the HLT menu the CPU scintillator module wires noiseMip as the
    # HGCAL_noise_heback PSet (scaleByDose etc.); the scalar MIP noise the SoA
    # producer needs lives in its noise_MIP field.
    def _scalar(x):
        if hasattr(x, "noise_MIP"):
            return x.noise_MIP.value()
        try:
            return x.value()
        except Exception:
            return None

    sci_noise = {}
    _cpuSci = getattr(process, "hltHgcalLayerClustersHSci", None)
    if _cpuSci is not None and hasattr(_cpuSci, "plugin"):
        _pl = _cpuSci.plugin
        if hasattr(_pl, "noiseMip"):
            _nm = _scalar(_pl.noiseMip)
            if _nm is not None:
                sci_noise["noiseMip"] = cms.double(_nm)
        if hasattr(_pl, "sciThicknessCorrection"):
            _sc = _scalar(_pl.sciThicknessCorrection)
            if _sc is not None:
                sci_noise["sciThicknessCorrection"] = cms.double(_sc)

    # FH and BH device chains, cloned from the EE modules.
    #   tag -> (detector, recHits, deltac)
    specs = [
        ("HSi", "FH", cms.InputTag("hltHGCalRecHit", "HGCHEFRecHits"), 1.3),
        ("HSci", "BH", cms.InputTag("hltHGCalRecHit", "HGCHEBRecHits"), 0.0315),
    ]
    new_modules = []
    for tag, det, recHits, deltac in specs:
        # BH needs the scintillator noise wired through so sigmaNoise != 0.
        rh_kwargs = dict(detector=det, recHits=recHits)
        if det == "BH":
            rh_kwargs.update(sci_noise)
        soarh = process.hltHgcalSoARecHitsProducer.clone(**rh_kwargs)
        # 'detector' is a new parameter on these modules, so pass it as cms.string.
        clue = process.hltHgcalSoARecHitsLayerClustersProducer.clone(
            hgcalRecHitsSoA="hltHgcalSoARecHits" + tag, deltac=deltac, detector=cms.string(det))
        agg = process.hltHgcalSoALayerClustersProducer.clone(
            hgcalRecHitsSoA="hltHgcalSoARecHits" + tag,
            hgcalRecHitsLayerClustersSoA="hltHgcalCLUEstering" + tag,
            detector=cms.string(det))
        conv = process.hltHgCalLayerClustersFromSoAProducer.clone(
            detector=det,
            src="hltHgcalSoALayerClusters" + tag,
            hgcalRecHitsSoA="hltHgcalSoARecHits" + tag,
            hgcalRecHitsLayerClustersSoA="hltHgcalCLUEstering" + tag)
        setattr(process, "hltHgcalSoARecHits" + tag, soarh)
        setattr(process, "hltHgcalCLUEstering" + tag, clue)
        setattr(process, "hltHgcalSoALayerClusters" + tag, agg)
        setattr(process, "hltHgCalLayerClustersFromSoAProducer" + tag, conv)
        new_modules += [soarh, clue, agg, conv]

    # Repoint the merge at the device converters for all three detectors.
    process.hltMergeLayerClusters.layerClusters = cms.VInputTag(
        "hltHgCalLayerClustersFromSoAProducer",       # EE
        "hltHgCalLayerClustersFromSoAProducerHSci",   # BH
        "hltHgCalLayerClustersFromSoAProducerHSi",    # FH
    )
    process.hltMergeLayerClusters.time_layerclusters = cms.VInputTag(
        "hltHgCalLayerClustersFromSoAProducer:timeLayerCluster",
        "hltHgCalLayerClustersFromSoAProducerHSci:timeLayerCluster",
        "hltHgCalLayerClustersFromSoAProducerHSi:timeLayerCluster",
    )

    # Schedule the new producers: associate a Task with every Path that runs the merge.
    process.hltHgcalCLUEsteringFHBHTask = cms.Task(*new_modules)
    for pn in process.paths_().keys():
        pth = getattr(process, pn)
        if "hltMergeLayerClusters" in pth.moduleNames():
            pth.associate(process.hltHgcalCLUEsteringFHBHTask)
    return process
