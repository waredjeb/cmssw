#include "RecoHGCal/TICL/interface/TrackstersLinkingPluginFactory.h"
#include "TrackstersLinkingbyPCA.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(TrackstersLinkingFactory, "TrackstersLinkingFactory");
DEFINE_EDM_VALIDATED_PLUGIN(TrackstersLinkingFactory, ticl::TrackstersLinkingbyPCA, "PCA");