#ifndef DataFormats_TICL_interface_alpaka_HitsAndFractionsDevice_h
#define DataFormats_TICL_interface_alpaka_HitsAndFractionsDevice_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/HitAndFraction.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  // Same layout as the host HitsAndFractionsHost (AssociationMap<int, HitAndFraction>),
  // so the device collection can be copied to/from the host collection directly.
  using HitsAndFractionsDevice = PortableCollection<::ticl::AssociationMap<int, ::ticl::HitAndFraction>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl

#endif
