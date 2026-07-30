import FWCore.ParameterSet.Config as cms

ceh_layerClusters = [
    "hltHgcalLayerClustersHSci",
    "hltHgcalLayerClustersHSi"
]

# Heterogeneous CE-H layer clusters, converted back from the SoA
ceh_layerClustersFromSoA = [
    "hltHgCalLayerClustersFromSoAProducerHSci",
    "hltHgCalLayerClustersFromSoAProducerHSi"
]

ceh_layerClustersFromSoASerialSync = [
    "hltHgCalLayerClustersFromSoAProducerHSciSerialSync",
    "hltHgCalLayerClustersFromSoAProducerHSiSerialSync"
]

barrel_layerClusters = [
    "hltBarrelLayerClustersEB",
    "hltBarrelLayerClustersHB"
]

# Define the producer with ceh lists
hltMergeLayerClusters = cms.EDProducer("MergeClusterProducer",
    layerClusters = cms.VInputTag("hltHgcalLayerClustersEE", *ceh_layerClusters),
)

hltMergeLayerClustersSerialSync = cms.EDProducer("MergeClusterProducer",
    layerClusters = cms.VInputTag("hltHgCalLayerClustersFromSoAProducerSerialSync", *ceh_layerClustersFromSoASerialSync),
)

# Process modifiers: ticl_barrel and alpaka
from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

(alpaka & ~ticl_barrel).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgCalLayerClustersFromSoAProducer", *ceh_layerClustersFromSoA]
)

(ticl_barrel & ~alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgcalLayerClustersEE", *ceh_layerClusters, *barrel_layerClusters]
)

(ticl_barrel & alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgCalLayerClustersFromSoAProducer", *ceh_layerClustersFromSoA, *barrel_layerClusters]
)
