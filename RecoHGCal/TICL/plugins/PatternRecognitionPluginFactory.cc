#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "PatternRecognitionbyCA.h"
#include "PatternRecognitionbyCLUE3D.h"
#include "PatternRecognitionbyFastJet.h"
#include "PatternRecognitionbyRecovery.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionFactory, "PatternRecognitionFactory");
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionHFNoseFactory, "PatternRecognitionHFNoseFactory");
// Note: PatternRecognitionFactoryPortable is registered in src/PatternRecognitionFactoryPortable.cc

DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyCA<ticl::TICLLayerTilesHost>, "CA");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory,
                            ticl::PatternRecognitionbyCLUE3D<ticl::TICLLayerTilesHost>,
                            "CLUE3D");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory,
                            ticl::PatternRecognitionbyFastJet<ticl::TICLLayerTilesHost>,
                            "FastJet");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory,
                            ticl::PatternRecognitionbyRecovery<ticl::TICLLayerTilesHost>,
                            "Recovery");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHFNoseFactory,
                            ticl::PatternRecognitionbyCA<ticl::TICLLayerTilesHFNoseHost>,
                            "CA");
