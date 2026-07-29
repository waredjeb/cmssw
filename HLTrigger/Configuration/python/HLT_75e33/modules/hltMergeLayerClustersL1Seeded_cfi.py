import FWCore.ParameterSet.Config as cms

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
