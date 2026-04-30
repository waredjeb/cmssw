#ifndef RecoHGCal_TICL_PatternRecognitionByCLUEsteringWrapper_H
#define RecoHGCal_TICL_PatternRecognitionByCLUEsteringWrapper_H

#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBasePortable.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUEstering.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <memory>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Wrapper class that adapts backend-specific implementation to backend-independent interface
  class PatternRecognitionByCLUEsteringWrapper final : public ::ticl::PatternRecognitionAlgoBasePortable {
  private:
    std::unique_ptr<PatternRecognitionAlgoBase> impl_;

  public:
    PatternRecognitionByCLUEsteringWrapper(const edm::ParameterSet& config)
        : ::ticl::PatternRecognitionAlgoBasePortable(config),
          impl_(std::make_unique<PatternRecognitionByCLUEstering>(config)) {}

    ~PatternRecognitionByCLUEsteringWrapper() override = default;

    // Provide access to backend-specific implementation
    PatternRecognitionAlgoBase* getImpl() { return impl_.get(); }

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      PatternRecognitionByCLUEstering::fillPSetDescription(iDesc);
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
