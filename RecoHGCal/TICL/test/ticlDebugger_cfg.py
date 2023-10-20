import FWCore.ParameterSet.Config as cms

process = cms.Process("TICLDEBUG")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T25', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data/wredjeb/Samples/step3/step3_16314654_39.root'
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

process.p = cms.Path(process.ticlDebugger+process.caloParticleDebugger)

