import FWCore.ParameterSet.Config as cms

hltHGCALGPUvsCPUComparisonHists = cms.EDProducer("HGCALGPUvsCPUComparisonHists",
                                                 monitoredLayerClusters = cms.InputTag("hltCaloClustersFromSoA"),
                                                 referenceLayerClusters = cms.InputTag("hltCaloClustersFromSoASerialSync"),
                                                 topFolderName = cms.string('HLT/HeterogeneousComparisons/HGCalMonitoring')
                                                 )
