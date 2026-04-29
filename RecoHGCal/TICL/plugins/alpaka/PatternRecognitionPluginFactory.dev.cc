#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "RecoHGCal/TICL/plugins/alpaka/PatternRecognitionByCLUEstering.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "FWCore/Utilities/interface/stringize.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <memory>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Wrapper class that adapts backend-specific implementation to backend-independent interface
  class PatternRecognitionByCLUEsteringWrapper final : public ticl::PatternRecognitionAlgoBasePortable {
  private:
    std::unique_ptr<PatternRecognitionAlgoBase> impl_;

  public:
    PatternRecognitionByCLUEsteringWrapper(const edm::ParameterSet& config)
        : ticl::PatternRecognitionAlgoBasePortable(config),
          impl_(std::make_unique<PatternRecognitionByCLUEstering>(config)) {}

    ~PatternRecognitionByCLUEsteringWrapper() override = default;

    // Provide access to backend-specific implementation
    PatternRecognitionAlgoBase* getImpl() { return impl_.get(); }

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      PatternRecognitionByCLUEstering::fillPSetDescription(iDesc);
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// Register with backend-specific name to avoid conflicts
// Serial backend registers as "alpaka_serial_sync::CLUEstering"
// CUDA backend registers as "alpaka_cuda_async::CLUEstering"
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactoryPortable,
                            ALPAKA_ACCELERATOR_NAMESPACE::PatternRecognitionByCLUEsteringWrapper,
                            EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) "::CLUEstering");
