import FWCore.ParameterSet.Config as cms

# The clusters themselves are already merged by hltHgcalSoALayerClustersProducer;
# what is rebuilt here is the hits-and-fractions map, which needs the
# per-subdetector rechits. The entries must therefore be listed in the same order
# as in that module, and clusterOffsets tells where each of them starts in the
# merged collection.
hltHgCalLayerClustersFromSoAProducer = cms.EDProducer("HGCalLayerClustersFromSoAProducer",
    src = cms.InputTag("hltHgcalSoALayerClustersProducer"),
    clusterOffsets = cms.InputTag("hltHgcalSoALayerClustersProducer", "clusterOffsets"),
    layerClusters = cms.VPSet(
        cms.PSet(
            detector = cms.string('EE'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducer"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducer")
        ),
        cms.PSet(
            detector = cms.string('BH'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSci"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSci")
        ),
        cms.PSet(
            detector = cms.string('FH'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSi"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSi")
        )
    ),
    nHitsTime = cms.uint32(3),
    timeClname = cms.string('timeLayerCluster')
)

hltHgCalLayerClustersFromSoAProducerSerialSync = cms.EDProducer("HGCalLayerClustersFromSoAProducer",
    src = cms.InputTag("hltHgcalSoALayerClustersProducerSerialSync"),
    clusterOffsets = cms.InputTag("hltHgcalSoALayerClustersProducerSerialSync", "clusterOffsets"),
    layerClusters = cms.VPSet(
        cms.PSet(
            detector = cms.string('EE'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerSerialSync")
        ),
        cms.PSet(
            detector = cms.string('BH'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSciSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSciSerialSync")
        ),
        cms.PSet(
            detector = cms.string('FH'),
            hgcalRecHitsLayerClustersSoA = cms.InputTag("hltHgcalSoARecHitsLayerClustersProducerHSiSerialSync"),
            hgcalRecHitsSoA = cms.InputTag("hltHgcalSoARecHitsProducerHSiSerialSync")
        )
    ),
    nHitsTime = cms.uint32(3),
    timeClname = cms.string('timeLayerCluster')
)
