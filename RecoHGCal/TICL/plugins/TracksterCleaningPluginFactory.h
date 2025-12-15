#ifndef RecoHGCal_TICL_TracksterCleaningPluginFactory_H
#define RecoHGCal_TICL_TracksterCleaningPluginFactory_H

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoHGCal/TICL/interface/TracksterCleaningAlgoBase.h"

typedef edmplugin::PluginFactory<ticl::TracksterCleaningAlgoBase*(const edm::ParameterSet&, edm::ConsumesCollector)>
    TracksterCleaningPluginFactory;

#endif