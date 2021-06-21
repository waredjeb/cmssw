#ifndef RecoHGCal_TICL_TrackstersLinkingPluginFactory_H
#define RecoHGCal_TICL_TrackstersLinkingPluginFactory_H

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoHGCal/TICL/interface/TrackstersLinkingAlgoBase.h"
#include "RecoHGCal/TICL/interface/GlobalCache.h"

typedef edmplugin::PluginFactory<ticl::TrackstersLinkingAlgoBaseT*(const edm::ParameterSet&, const ticl::CacheBase*)>
    TrackstersLinkingFactory;
#endif