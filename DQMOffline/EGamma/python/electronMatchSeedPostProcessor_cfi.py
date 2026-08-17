import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

electronMatchSeedPostProcessor = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring("HLT/EGamma/PixelSeedMatching"),
    efficiency = cms.vstring(
        #reconstruction efficiency: numerator and denominator both gen-binned,
        #so a matched electron lands in the same bin in both
        "'Eff_vs_pt' 'pixel seed matching efficiency vs p_{T};p_{T} [GeV];efficiency' eG_pt electronSim",
        "'Eff_vs_eta' 'pixel seed matching efficiency vs #eta;#eta;efficiency' eG_eta electronSim_eta",
        "'Eff_vs_phi' 'pixel seed matching efficiency vs #phi;#phi;efficiency' eG_phi electronSim_phi",
        #matched electrons: both sides from the same matched
        #population in the same (reco) variable
        #the three categories are exclusive and cover every seed size, so they sum to 1.
        #eQuad is only populated when the seeds carry >=4 hits, e.g. with the
        #hltEgammaPixelTrackSeeding process modifier
        "'FracFromDoublet_vs_pt' 'fraction of matched electrons from doublet seeds;p_{T} [GeV];fraction' eDouble eMatched",
        "'FracFromTriplet_vs_pt' 'fraction of matched electrons from triplet seeds;p_{T} [GeV];fraction' eTrip eMatched",
        "'FracFromQuadPlus_vs_pt' 'fraction of matched electrons from #geq4-hit seeds;p_{T} [GeV];fraction' eQuad eMatched",
    ),
    resolution = cms.vstring(),
    efficiencyProfile = cms.untracked.vstring(),
    verbose = cms.untracked.uint32(2),
    outputFileName = cms.untracked.string(""),
)