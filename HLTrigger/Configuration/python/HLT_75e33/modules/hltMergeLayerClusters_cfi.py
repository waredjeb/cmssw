import FWCore.ParameterSet.Config as cms

# HGCal layer clusters come from the device chain: the SoA->legacy converters
# emit the same legacy products the CPU producers used to, so the merge reads
# them directly. The clustering backend (CUDA, ROCm, serial) is chosen by the
# alpaka services, so no ProcessModifier is involved here.
hgcal_layerClusters = [
    "hltHgCalLayerClustersFromSoAProducer",
    "hltHgCalLayerClustersFromSoAProducerHSci",
    "hltHgCalLayerClustersFromSoAProducerHSi"
]
hgcal_time_layerClusters = [x + ":timeLayerCluster" for x in hgcal_layerClusters]

hgcal_layerClustersSerialSync = [x + "SerialSync" for x in hgcal_layerClusters]
hgcal_time_layerClustersSerialSync = [x + "SerialSync:timeLayerCluster" for x in hgcal_layerClusters]

barrel_layerClusters = [
    "hltBarrelLayerClustersEB",
    "hltBarrelLayerClustersHB"
]
barrel_time_layerClusters = [x + ":timeLayerCluster" for x in barrel_layerClusters]

hltMergeLayerClusters = cms.EDProducer("MergeClusterProducer",
    layerClusters = cms.VInputTag(*hgcal_layerClusters),
    time_layerclusters = cms.VInputTag(*hgcal_time_layerClusters),
)

hltMergeLayerClustersSerialSync = cms.EDProducer("MergeClusterProducer",
    layerClusters = cms.VInputTag(*hgcal_layerClustersSerialSync),
    time_layerclusters = cms.VInputTag(*hgcal_time_layerClustersSerialSync),
)

from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

ticl_barrel.toModify(hltMergeLayerClusters,
    layerClusters = [*hgcal_layerClusters, *barrel_layerClusters],
    time_layerclusters = [*hgcal_time_layerClusters, *barrel_time_layerClusters]
)
