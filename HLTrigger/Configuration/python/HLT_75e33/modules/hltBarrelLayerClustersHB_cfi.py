import FWCore.ParameterSet.Config as cms

hltBarrelLayerClustersHB = cms.EDProducer("BarrelLayerClustersFromSoAProducer",
    src = cms.InputTag("hltBarrelSoALayerClustersHB"),
    recHits = cms.InputTag("hltBarrelPFRecHitsSoAHB"),
    barrelRecHitsLayerClustersSoA = cms.InputTag("hltBarrelCLUEsteringHB"),
    detector = cms.string("HB"),
    nHitsTime = cms.uint32(3),
    timeClname = cms.string("timeLayerCluster"),
    timeResolutionCalc = cms.PSet(
        noiseTerm = cms.double(1.10889),
        constantTerm = cms.double(0.428192),
        corrTermLowE = cms.double(0.0510871),
        threshLowE = cms.double(0.5),
        constantTermLowE = cms.double(0.0),
        noiseTermLowE = cms.double(1.31883),
        threshHighE = cms.double(5.0)
    )
)
