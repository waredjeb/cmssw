#ifndef RecoHGCal_TICL_PatternRecognitionAlgoBasePortable_H
#define RecoHGCal_TICL_PatternRecognitionAlgoBasePortable_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Backend-independent base class (outside ALPAKA_ACCELERATOR_NAMESPACE)
// This allows plugin registration to happen once per backend with unique names
namespace ticl {

  class PatternRecognitionAlgoBasePortable {
  public:
    PatternRecognitionAlgoBasePortable(const edm::ParameterSet& conf) {}
    virtual ~PatternRecognitionAlgoBasePortable() = default;

    // Pure virtual interface - backend-specific implementations will override
    // Note: cannot use Queue or device-specific types here
    // Derived classes will provide backend-specific interfaces
  };

}  // namespace ticl

#endif
