#ifndef DataFormats_TICL_interface_HitAndFraction_h
#define DataFormats_TICL_interface_HitAndFraction_h

#include "DataFormats/DetId/interface/DetId.h"

namespace ticl {

  struct HitAndFraction {
    DetId hit;
    float fraction;
  };

}  // namespace ticl

#endif
