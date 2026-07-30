import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

# Scintillator tiles cluster in (eta, phi), so the critical distance needs the
# scintillator scale rather than the silicon (cm) one: deltac = 0.0315 and
# outlierDeltaFactor = 2.0 reproduce the CPU SciCLUE deltac/deltas/deltao of
# 0.0315/0.0315/0.063.
hltHgcalSoARecHitsLayerClustersProducerHSci = cms.EDProducer("HGCalSoARecHitsLayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    detector = cms.string('BH'),
    hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSci"),
    deltac = cms.double(0.0315),
    kappa = cms.double(9),
    outlierDeltaFactor = cms.double(2.0)
)

hltHgcalSoARecHitsLayerClustersProducerHSciSerialSync = makeSerialClone(hltHgcalSoARecHitsLayerClustersProducerHSci,
                                                                        hgcalRecHitsSoA = "hltHgcalSoARecHitsProducerHSciSerialSync"
)
