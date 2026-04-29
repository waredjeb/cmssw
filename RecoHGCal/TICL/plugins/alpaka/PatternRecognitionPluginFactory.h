#ifndef RecoHGCal_TICL_PatternRecognitionPluginFactory_Alpaka_H
#define RecoHGCal_TICL_PatternRecognitionPluginFactory_Alpaka_H

// This file is kept for backward compatibility but is no longer used.
// Pattern recognition algorithms now register with the backend-independent
// PatternRecognitionFactoryPortable defined in RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h
//
// Each Alpaka backend registers its implementations with backend-specific names
// (e.g., "alpaka_serial_sync::CLUEstering", "alpaka_cuda_async::CLUEstering")
// to avoid plugin name conflicts.

#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"

#endif
