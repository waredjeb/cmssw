#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUEsteringWrapper.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "FWCore/Utilities/interface/stringize.h"

// Register with backend-specific name to avoid conflicts
// Serial backend registers as "alpaka_serial_sync::CLUEstering"
// CUDA backend registers as "alpaka_cuda_async::CLUEstering"
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactoryPortable,
                            ALPAKA_ACCELERATOR_NAMESPACE::PatternRecognitionByCLUEsteringWrapper,
                            EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) "::CLUEstering");
