#ifndef DataFormats_TICL_interface_alpaka_HitsAndFractionsDevice_h
#define DataFormats_TICL_interface_alpaka_HitsAndFractionsDevice_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/HitAndFraction.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl {

  using HitsAndFractionsDevice = PortableCollection<::ticl::AssociationMap<int, ::ticl::HitAndFraction>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ticl

#endif
