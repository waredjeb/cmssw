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

        #'file:/afs/cern.ch/work/w/wredjeb/public/TracksterLinking/CMSSW_12_0_X_2021-06-20-2300/src/37888.0_SinglePiPt25Eta1p7_2p7+2026D84+SinglePiPt25Eta1p7_2p7_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_pi0_300.root'
        'file:/afs/cern.ch/work/w/wredjeb/public/TracksterLinking/CMSSW_12_0_X_2021-06-29-2300/src/37888.0_SinglePiPt25Eta1p7_2p7+2026D84+SinglePiPt25Eta1p7_2p7_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_pi0_300.root'
    )
)

# process.load("RecoHGCal.TICL.tracksterLinkingDebugger_cfi")
process.load("SimGeneral.Debugging.caloParticleDebugger_cfi")

process.tracksterLinkingDebugger = cms.EDAnalyzer('TrackstersLinkingDebugger',
  tracks = cms.InputTag('generalTracks'),
  tracksters = cms.InputTag('ticlTrackstersMerge'),
  caloparticles = cms.InputTag('mix', 'MergedCaloTruth'),
  min_angle = cms.double(0.0),
  do2Ddot = cms.bool(False),
  mightGet = cms.optional.untracked.vstring
)

# MessageLogger customizations
process.MessageLogger.cerr.enable = False
process.MessageLogger.cout.enable = False
label = 'TrackstersLinkingDebugger'
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

process.p = cms.Path(process.tracksterLinkingDebugger+process.caloParticleDebugger)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("output_pi0_gen_cos0.root"),
      closeFileFast = cms.untracked.bool(True)
  )


