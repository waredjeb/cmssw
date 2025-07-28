import FWCore.ParameterSet.Config as cms
from PhysicsTools.PyTorch.options import args
from PhysicsTools.PyTorch.modules import (
    torchtest_DataProducer_alpaka,
    torchtest_JitClassificationProducer_alpaka,
    torchtest_JitRegressionProducer_alpaka,
    torchtest_AotClassificationProducer_alpaka,
    torchtest_AotRegressionProducer_alpaka,
    torchtest_CombinatoricsProducer_alpaka,
)
  
args.parseArguments()
process = cms.Process("TestPyTorchHeterogeneousPipeline")

# enable multithreading
process.options.numberOfThreads = args.numberOfThreads if args.numberOfThreads > 1 else 1 
process.options.numberOfStreams = args.numberOfStreams if args.numberOfStreams > 1 else 1 

# enable alpaka and GPU support
process.load("Configuration.StandardSequences.Accelerators_cff")

# process a limited number of events
process.maxEvents.input = args.numberOfEvents if args.numberOfEvents > 1 else 1 

# empty source
process.source = cms.Source("EmptySource")

# print a message every event
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# do not print the time and trigger reports at the end of the job
process.options.wantSummary = False

# setup chain configs
process.DataProducer = torchtest_DataProducer_alpaka(
    batchSize = cms.uint32(args.batchSize if args.batchSize > 1 else 1),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)

# JIT models
process.JitClassificationProducer = torchtest_JitClassificationProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(args.classificationModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)
process.JitClassificationProducerCpu = torchtest_JitClassificationProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(args.classificationModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("serial_sync")  # force serial backend to emulate heterogeneous execution
    ),
)
process.JitRegressionProducer = torchtest_JitRegressionProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(args.regressionModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)
process.JitRegressionProducerCpu = torchtest_JitRegressionProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(args.regressionModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string("serial_sync")  # force serial backend to emulate heterogeneous execution
    ),
)

# AOT models
regressionModelPath = args.regressionModelPathCpu
if args.backend == "cuda_async":
    regressionModelPath = args.regressionModelPathCuda
process.AotRegressionProducer = torchtest_AotRegressionProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(regressionModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)
classificationModelPath = args.classificationModelPathCpu
if args.backend == "cuda_async":
    classificationModelPath = args.classificationModelPathCuda
process.AotClassificationProducer = torchtest_AotClassificationProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    modelPath = cms.FileInPath(classificationModelPath),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)

process.CombinatoricsProducer = torchtest_CombinatoricsProducer_alpaka(
    inputs = cms.InputTag('DataProducer'),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
)

# schedule the modules
process.path = cms.Path(
    process.DataProducer +
    process.JitClassificationProducer +
    # process.JitClassificationProducerCpu +
    # process.JitRegressionProducer +
    process.JitRegressionProducerCpu +

    # disable since no proper support is available yet
    # process.AotRegressionProducer +
    # process.AotClassificationProducer +

    process.CombinatoricsProducer
)
