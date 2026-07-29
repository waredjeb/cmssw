import FWCore.ParameterSet.Config as cms

# Layer-cluster timing now travels inside the cluster SoA, so the merge no
# longer takes a parallel list of time ValueMaps.
hltMergeLayerClustersL1Seeded = cms.EDProducer("MergeClusterProducer",
    layerClusters = cms.VInputTag("hltHgcalLayerClustersEEL1Seeded", "hltHgcalLayerClustersHSciL1Seeded", "hltHgcalLayerClustersHSiL1Seeded"),
    mightGet = cms.optional.untracked.vstring
)

from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

layerClusters = ["hltHgcalLayerClustersEEL1Seeded",
                 "hltHgcalLayerClustersHSciL1Seeded",
                 "hltHgcalLayerClustersHSiL1Seeded",
                 "hltBarrelLayerClustersEBL1Seeded"]

ticl_barrel.toModify(hltMergeLayerClustersL1Seeded, layerClusters = layerClusters)
