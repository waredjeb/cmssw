#ifndef DataFormats_HGCalReco_interface_TrackstersSoA_h
#define DataFormats_HGCalReco_interface_TrackstersSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

GENERATE_SOA_LAYOUT(TrackstersSoALayout,
                    // columns: one value per element
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z),
                    SOA_COLUMN(float, energy)
)

using TrackstersSoA = TrackstersSoALayout<>;

#endif
