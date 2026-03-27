
#pragma once

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/TICL/interface/HitAndFraction.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace ticl {

  // Simple SoA layout for hits and fractions - no custom view methods
  // This avoids the PortableHostCollection + complex BLOCKS incompatibility
  GENERATE_SOA_LAYOUT(HitsAndFractionsSoALayout,
    SOA_COLUMN(int, offsets),           // Offset for each cluster
    SOA_COLUMN(HitAndFraction, values)  // All hit-fraction pairs
  )

  using HitsAndFractionsSoA = HitsAndFractionsSoALayout<>;
  using HitsAndFractionsHost = PortableHostCollection<HitsAndFractionsSoA>;

}  // namespace ticl
