import FWCore.ParameterSet.Config as cms

process = cms.Process("TICLDEBUG")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:step3.root'
        # 'file:/afs/cern.ch/work/w/wredjeb/public/TICLv4/Associators/CMSSW_12_4_0_pre2/src/38694.203_CloseByPGun_CE_E_Front_300um+2026D86_ticl_v4+CE_E_Front_300um_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_4cm.root'
        # 'file:/afs/cern.ch/work/w/wredjeb/public/TICLv4/Associators/CMSSW_12_4_0_pre2/src/38694.203_CloseByPGun_CE_E_Front_300um+2026D86_ticl_v4+CE_E_Front_300um_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_4cm_sym.root'
        'file:/afs/cern.ch/work/w/wredjeb/public/TICLv4/Associators/CMSSW_12_4_0_pre2/src/38694.203_CloseByPGun_CE_E_Front_300um+2026D86_ticl_v4+CE_E_Front_300um_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_1p5cm_asym.root'
        # 'file:/data2/user/asavona/CMSSW_12_3_0_pre5/src/SampleProduction/38694.201_CloseByPGun_AsymEnergy_HD/38694.201_CloseByPGun_CE_E_HD_AsymEnergy_Delta4_0_Front_300um+2026D86_ticl_clue3D+CE_E_Front_300um_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root'
    )
)

process.load("RecoHGCal.TICL.ticlDebugger_cfi")
process.load("SimGeneral.Debugging.caloParticleDebugger_cfi")


# MessageLogger customizations
process.MessageLogger.cerr.enable = False
process.MessageLogger.cout.enable = False
label = 'TICLDebugger'
messageLogger = dict()
main_key = '%sMessageLogger'%(label)
messageLogger[main_key] = dict(
        filename = '%s.log' % (label),
        threshold = 'INFO',
        default = dict(limit=0)
        )
messageLogger[main_key][label] = dict(limit=-1)
# First create defaults
setattr(process.MessageLogger.files, label, dict())
# Then modify them
setattr(process.MessageLogger.files, label, messageLogger[main_key])

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo_1p5cm_asym.root") )

process.p = cms.Path(process.ticlDebugger+process.caloParticleDebugger)

