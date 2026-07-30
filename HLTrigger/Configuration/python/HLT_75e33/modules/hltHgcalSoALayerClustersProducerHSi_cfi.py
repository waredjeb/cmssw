import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

hltHgcalSoALayerClustersProducerHSi = cms.EDProducer("HGCalSoALayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    detector = cms.string('FH'),
    hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSi", "maxLayerPerSide"),
    hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSi"),
    hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSi"),
    positionDeltaRho2 = cms.double(1.69),
    thresholdW0 = cms.double(2.9)
)

hltHgcalSoALayerClustersProducerHSiSerialSync = makeSerialClone(hltHgcalSoALayerClustersProducerHSi,
                                                                #feed the upstream serial modules in
                                                                hgcalMaxLayerPerSide = ("hltHgcalSoARecHitsProducerHSiSerialSync", "maxLayerPerSide"),
                                                                hgcalRecHitsLayerClustersSoA = "hltHgcalSoARecHitsLayerClustersProducerHSiSerialSync",
                                                                hgcalRecHitsSoA = "hltHgcalSoARecHitsProducerHSiSerialSync"
)
