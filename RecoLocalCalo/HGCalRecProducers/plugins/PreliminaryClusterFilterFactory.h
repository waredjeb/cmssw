#ifndef RecoHGCal_TICL_PreliminaryClusterFilterFactory_H
#define RecoHGCal_TICL_PreliminaryClusterFilterFactory_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "PreliminaryClusterFilterBase.h"
#include "PreliminaryClusterFilterByCP.h"

typedef edmplugin::PluginFactory<ticl::PreliminaryClusterFilterBase*(const edm::ParameterSet&)> PreliminaryClusterFilterFactory;

#endif
