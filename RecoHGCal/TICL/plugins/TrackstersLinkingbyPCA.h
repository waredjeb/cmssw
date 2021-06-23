#ifndef __RecoHGCal_TICL_PRbyCA_H__
#define __RecoHGCal_TICL_PRbyCA_H__
#include <memory>  // unique_ptr
#include "RecoHGCal/TICL/interface/TrackstersLinkingAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

namespace ticl {
  class TrackstersLinkingbyPCA final : public TrackstersLinkingAlgoBaseT {
  public:
    TrackstersLinkingbyPCA(const edm::ParameterSet& conf, const CacheBase* cache);
    ~TrackstersLinkingbyPCA() override;

    void trackstersInfo(const typename TrackstersLinkingAlgoBaseT::Inputs& input);
    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {}
  };
}  // namespace ticl
#endif