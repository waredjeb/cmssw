#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "TracksterCleaningByBeta.h"
#include "RecoHGCal/TICL/plugins/TracksterCleaningPluginFactory.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(TracksterCleaningPluginFactory, "TracksterCleaningPluginFactory");

DEFINE_EDM_VALIDATED_PLUGIN(TracksterCleaningPluginFactory, 
                            ticl::TracksterCleaningByBeta, 
                            "Beta");