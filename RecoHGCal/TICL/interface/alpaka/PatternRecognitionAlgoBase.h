
#pragma once

#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  namespace ticl = ::ticl;

  class PatternRecognitionAlgoBase {
  protected:
    int m_verbosity;

  public:
    PatternRecognitionAlgoBase(const edm::ParameterSet& conf) : m_verbosity(conf.getParameter<int>("verbosity")) {}
    virtual ~PatternRecognitionAlgoBase() = default;

    virtual void makeTracksters(Queue& queue,
                                const HGCalSoAClustersDeviceCollection& layerClusters,
                                std::vector<ticl::Trackster>& result) = 0;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl
