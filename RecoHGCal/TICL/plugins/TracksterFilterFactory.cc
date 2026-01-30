// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 01/2026

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "TracksterFilterFactory.h"

#include "TracksterFilterByPDGID.h"

using namespace ticl;

EDM_REGISTER_PLUGINFACTORY(TracksterFilterFactory, "TracksterFilterFactory");
DEFINE_EDM_PLUGIN(TracksterFilterFactory, TracksterFilterByPDGID, "TracksterFilterByPDGID");
