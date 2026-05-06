#ifndef RecoHGCal_TICL_PatternRecognitionByCLUE3DWrapper_H
#define RecoHGCal_TICL_PatternRecognitionByCLUE3DWrapper_H

#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBasePortable.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUE3D.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <memory>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Wrapper class that adapts backend-specific implementation to backend-independent interface
  class PatternRecognitionByCLUE3DWrapper final : public ::ticl::PatternRecognitionAlgoBasePortable {
  private:
    std::unique_ptr<PatternRecognitionAlgoBase> impl_;

  public:
    PatternRecognitionByCLUE3DWrapper(const edm::ParameterSet& config)
        : ::ticl::PatternRecognitionAlgoBasePortable(config),
          impl_(std::make_unique<PatternRecognitionByCLUE3D>(config)) {}

    ~PatternRecognitionByCLUE3DWrapper() override = default;

    // Provide access to backend-specific implementation
    PatternRecognitionAlgoBase* getImpl() { return impl_.get(); }

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      PatternRecognitionByCLUE3D::fillPSetDescription(iDesc);
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
