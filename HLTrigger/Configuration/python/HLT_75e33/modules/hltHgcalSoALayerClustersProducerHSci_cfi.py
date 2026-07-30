import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

# thresholdW0 and positionDeltaRho2 only enter the silicon position recipe;
# the scintillator (BH) position is a plain energy-weighted mean of all cells.
hltHgcalSoALayerClustersProducerHSci = cms.EDProducer("HGCalSoALayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    detector = cms.string('BH'),
    hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSci", "maxLayerPerSide"),
    hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSci"),
    hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSci"),
    positionDeltaRho2 = cms.double(1.69),
    thresholdW0 = cms.double(2.9)
)

hltHgcalSoALayerClustersProducerHSciSerialSync = makeSerialClone(hltHgcalSoALayerClustersProducerHSci,
                                                                 #feed the upstream serial modules in
                                                                 hgcalMaxLayerPerSide = ("hltHgcalSoARecHitsProducerHSciSerialSync", "maxLayerPerSide"),
                                                                 hgcalRecHitsLayerClustersSoA = "hltHgcalSoARecHitsLayerClustersProducerHSciSerialSync",
                                                                 hgcalRecHitsSoA = "hltHgcalSoARecHitsProducerHSciSerialSync"
)
