import FWCore.ParameterSet.Config as cms

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

# Process modifiers: ticl_barrel and alpaka
from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

# With alpaka, hltHgCalLayerClustersFromSoAProducer already emits EE, HSci and HSi
# as a single collection: nothing is left to merge unless the barrel is also run,
# so the module is only kept for ticl_barrel and dropped from the other sequences.
(ticl_barrel & ~alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgcalLayerClustersEE", *ceh_layerClusters, *barrel_layerClusters]
)

(ticl_barrel & alpaka).toModify(hltMergeLayerClusters,
    layerClusters = ["hltHgCalLayerClustersFromSoAProducer", *barrel_layerClusters]
)
