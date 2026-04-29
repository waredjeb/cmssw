// Factory registration for backend-independent Alpaka pattern recognition
// This file is in the main TICL library so the factory is available to all plugins

#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"

// Register backend-independent Alpaka factory ONCE in the main library
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionFactoryPortable, "PatternRecognitionFactoryPortable");
