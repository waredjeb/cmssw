#include <algorithm>
#include <set>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TrackstersLinkingbyPCA.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"

using namespace ticl;
TrackstersLinkingbyPCA::TrackstersLinkingbyPCA(const edm::ParameterSet &conf, const CacheBase *cache)
    : TrackstersLinkingAlgoBaseT(conf, cache){};

TrackstersLinkingbyPCA::~TrackstersLinkingbyPCA(){};