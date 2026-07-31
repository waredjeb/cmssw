import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

# The clustering runs per subdetector, but the clusters of all of them are
# assembled into a single SoA, each subdetector filling a slice of it through
# its own view. The order of the entries below is the order of the clusters in
# that collection, and is kept identical to the one hltMergeLayerClusters used
# (EE, HSci, HSi), so that the layer cluster indices are unchanged.
#
hltHgcalSoALayerClustersProducer = cms.EDProducer("HGCalSoALayerClustersProducer@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    layerClusters = cms.VPSet(
        cms.PSet(
            detector = cms.string('EE'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducer", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducer"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducer"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        ),
        cms.PSet(
            detector = cms.string('BH'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSci", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSci"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSci"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        ),
        cms.PSet(
            detector = cms.string('FH'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSi", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSi"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSi"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        )
    )
)

hltHgcalSoALayerClustersProducerSerialSync = makeSerialClone(hltHgcalSoALayerClustersProducer,
    #feed the upstream serial modules in
    layerClusters = cms.VPSet(
        cms.PSet(
            detector = cms.string('EE'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerSerialSync", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerSerialSync"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        ),
        cms.PSet(
            detector = cms.string('BH'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSciSerialSync", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSciSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSciSerialSync"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        ),
        cms.PSet(
            detector = cms.string('FH'),
            hgcalMaxLayerPerSide = cms.InputTag("hltHgcalSoARecHitsProducerHSiSerialSync", "maxLayerPerSide"),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSiSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSiSerialSync"),
            positionDeltaRho2 = cms.double(1.69),
            thresholdW0 = cms.double(2.9)
        )
    )
)
