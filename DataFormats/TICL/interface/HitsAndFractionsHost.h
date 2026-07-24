#ifndef DataFormats_TICL_interface_HitsAndFractionsHost_h
#define DataFormats_TICL_interface_HitsAndFractionsHost_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/HitAndFraction.h"

namespace ticl {

  // Association map keyed by layer-cluster index, valued by the list of
  // {DetId, fraction} that belong to each cluster (the "hits and fractions").
  // AssociationMap<int, HitAndFraction> resolves to Layout<128, false>, which is
  // identical to the default CacheLineSize::defaultSize (128) / relaxed layout,
  // so host and device collections share the exact same layout type.
  using TICLAssociationMap_t = AssociationMap<int, HitAndFraction>;
  using HitsAndFractionsHost = PortableHostCollection<TICLAssociationMap_t>;

}  // namespace ticl

#endif
