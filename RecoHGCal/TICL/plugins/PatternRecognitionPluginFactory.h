#ifndef RecoHGCal_TICL_PatternRecognitionPluginFactory_H
#define RecoHGCal_TICL_PatternRecognitionPluginFactory_H

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoHGCal/TICL/interface/GlobalCache.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBasePortable.h"

typedef edmplugin::PluginFactory<ticl::PatternRecognitionAlgoBaseT<ticl::TICLLayerTilesHost>*(const edm::ParameterSet&,
                                                                                              edm::ConsumesCollector)>
    PatternRecognitionFactory;
typedef edmplugin::PluginFactory<ticl::PatternRecognitionAlgoBaseT<ticl::TICLLayerTilesHFNoseHost>*(
    const edm::ParameterSet&, edm::ConsumesCollector)>
    PatternRecognitionHFNoseFactory;

// Backend-independent factory for Alpaka pattern recognition
// Each backend registers with a unique name (including backend namespace)
typedef edmplugin::PluginFactory<ticl::PatternRecognitionAlgoBasePortable*(const edm::ParameterSet&)>
    PatternRecognitionFactoryPortable;

#endif
