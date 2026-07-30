import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

# CE-H silicon uses the same CLUE parameters as CE-E (distances in cm); they
# are kept as a separate instance so that they can be tuned independently.
hltHgcalSoARecHitsLayerClustersProducerHSi = cms.EDProducer("HGCalSoARecHitsLayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    detector = cms.string('FH'),
    hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSi"),
    deltac = cms.double(1.3),
    kappa = cms.double(9),
    outlierDeltaFactor = cms.double(2.0)
)

hltHgcalSoARecHitsLayerClustersProducerHSiSerialSync = makeSerialClone(hltHgcalSoARecHitsLayerClustersProducerHSi,
                                                                       hgcalRecHitsSoA = "hltHgcalSoARecHitsProducerHSiSerialSync"
)
