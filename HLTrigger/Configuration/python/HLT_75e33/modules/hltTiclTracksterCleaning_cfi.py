import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.tracksterCleaningProducer_cfi import tracksterCleaningProducer as _tracksterCleaningProducer

popt1 = dict(
    betaContamMin = 1.3649,
    R0            = 0.30,
    sigmaZ        = 12.28,
    sigmaT        = 0.0663,
    sigmaDR       = 0.02,
    zPower        = 1.0,
    tPower        = 2.0,
    drPower       = 1.0,     
    wmin          = 0.0001,
)

popt2 = dict(
    betaContamMin = 1.259,
    R0            = 0.25,
    sigmaZ        = 5.,
    sigmaT        = 0.03,
    sigmaDR       = 0.0941,
    zPower        = 2.0,
    tPower        = 0.0,
    drPower       = 0.0,     
    wmin          = 0.001,
)

popt3 = dict(
    betaContamMin = 2.48,
    R0            = 0.1,
    sigmaZ        = 5.,
    sigmaT        = 0.09,
    sigmaDR       = 0.2,
    zPower        = 2.0,
    tPower        = 0.0,
    drPower       = 1.0,     
    wmin          = 0.0007,
)

hltTiclTracksterCleaning = _tracksterCleaningProducer.clone(
    linkedTracksters      = "hltTiclTracksterLinks",
    clue3DTracksters      = "hltTiclTrackstersCLUE3DHigh",
    layer_clusters        = "hltMergeLayerClusters",
    clue3DInLinkedIndices = ("hltTiclTracksterLinks", "linkedTracksterIdToInputTracksterId"),
    algo_verbosity = 0,
    cleaner = dict(
        type = "Beta",
        algo_verbosity = 1,
        useRawEnergy = True,
        epsE = 1e-6,
        epsDR = 1e-6,
        weightMode = True,
        emitDroppedAsStandalone = False,
        zAbsCut = 25,
        tAbsCut = 0.15,
        doPruning = False,

        betaContamMin = popt3["betaContamMin"],
        R0            = popt3["R0"],
        sigmaZ        = popt3["sigmaZ"],
        sigmaT        = popt3["sigmaT"],
        sigmaDR       = popt3["sigmaDR"],
        zPower        = popt3["zPower"],
        tPower        = popt3["tPower"],
        drPower       = popt3["drPower"],
        wmin          = popt3["wmin"],

        pruneWmin = 0.01,
        pruneUseSeparateKernels = False,
        sigmaZ_prune = 12.5,
        sigmaT_prune = 0.08,
        sigmaDR_prune = 0.08,
        zPower_prune = 1.5,
        tPower_prune = 0.5,
        drPower_prune = 0.5,
    ),
)