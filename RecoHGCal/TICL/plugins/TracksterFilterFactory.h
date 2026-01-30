// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 01/2026

#ifndef RecoHGCal_TICL_TracksterFilterFactory_h
#define RecoHGCal_TICL_TracksterFilterFactory_h

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "RecoHGCal/TICL/plugins/TracksterFilterBase.h"

using TracksterFilterFactory = edmplugin::PluginFactory<ticl::TracksterFilterBase*(const edm::ParameterSet&)>;

#endif
