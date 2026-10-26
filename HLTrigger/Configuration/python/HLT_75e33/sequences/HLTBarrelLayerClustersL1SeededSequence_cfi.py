import FWCore.ParameterSet.Config as cms

# L1-seeded EB layer clustering (CLUEstering), on the ECAL rechits in the L1 e/gamma regions.
# Same chain as in HLTBarrelLayerClustersSequence.
from ..modules.hltBarrelRecHitsSoAEBL1Seeded_cfi import *
from ..modules.hltBarrelPFRecHitsSoAEBL1Seeded_cfi import *
from ..modules.hltBarrelCLUEsteringEBL1Seeded_cfi import *
from ..modules.hltBarrelSoALayerClustersEBL1Seeded_cfi import *
from ..modules.hltBarrelLayerClustersEBL1Seeded_cfi import *

HLTBarrelLayerClustersL1SeededSequence = cms.Sequence(
        hltBarrelRecHitsSoAEBL1Seeded+
        hltBarrelPFRecHitsSoAEBL1Seeded+
        hltBarrelCLUEsteringEBL1Seeded+
        hltBarrelSoALayerClustersEBL1Seeded+
        hltBarrelLayerClustersEBL1Seeded)
