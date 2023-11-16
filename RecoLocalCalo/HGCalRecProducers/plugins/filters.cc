#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "PreliminaryClusterFilterFactory.h"

#include "PreliminaryClusterFilterByCP.h"

using namespace ticl;

DEFINE_EDM_PLUGIN(PreliminaryClusterFilterFactory, PreliminaryClusterFilterByCP, "PreliminaryClusterFilterByCP");
