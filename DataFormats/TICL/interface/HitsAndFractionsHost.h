#ifndef DataFormats_TICL_interface_HitsAndFractionsHost_h
#define DataFormats_TICL_interface_HitsAndFractionsHost_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/HitAndFraction.h"

namespace ticl {

  using TICLAssociationMap_t = AssociationMapLayout<int, HitAndFraction>::Layout<128, false>;
  using HitsAndFractionsHost = PortableHostCollection<TICLAssociationMap_t>;

}  // namespace ticl

#endif
