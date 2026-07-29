import FWCore.ParameterSet.Config as cms

# Layer-cluster timing now travels inside the cluster SoA, so the merge no
# longer takes a parallel list of time ValueMaps.

ceh_layerClusters = [
    "hltHgcalLayerClustersHSci",
    "hltHgcalLayerClustersHSi"
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
    layerClusters = cms.VInputTag("hltHgCalLayerClustersFromSoAProducerSerialSync", *ceh_layerClusters),
)

# Process modifiers: ticl_barrel and alpaka
from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

(alpaka & ~ticl_barrel).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgCalLayerClustersFromSoAProducer", *ceh_layerClusters]
)

(ticl_barrel & ~alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgcalLayerClustersEE", *ceh_layerClusters, *barrel_layerClusters]
)

(ticl_barrel & alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgCalLayerClustersFromSoAProducer", *ceh_layerClusters, *barrel_layerClusters]
)
