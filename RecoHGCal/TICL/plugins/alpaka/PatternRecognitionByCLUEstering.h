
#pragma once

#include "CondCore/CondDB/interface/Exception.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <algorithm>
#include <array>
#include <ranges>
#include <unordered_map>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PatternRecognitionByCLUEstering final : public PatternRecognitionAlgoBase {
  private:
    double m_rhoc;
    double m_dc;
    double m_dm;

  public:
    PatternRecognitionByCLUEstering(const edm::ParameterSet& config)
        : PatternRecognitionAlgoBase(config), m_rhoc(config.getParameter<double>("rho_c")) {
      auto dc = config.getParameter<double>("dc");
      auto dm = config.getParameter<double>("dm");
      m_dc = dc; 
      m_dm = dm;
    }
    ~PatternRecognitionByCLUEstering() override = default;

    void makeTracksters(Queue& queue,
                        const HGCalSoAClustersDeviceCollection& lc,
                        std::vector<ticl::Trackster>& tracksters) override;

    static void fillPSetDescription(::edm::ParameterSetDescription& iDesc);
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
