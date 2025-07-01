#ifndef DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h
#define DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/HGCalReco/interface/TrackstersSoA.h"

using TrackstersSoAHostCollection = PortableHostCollection<GNNNodeSoA>;
using TrackstersEdgeSoAHostCollection = PortableHostCollection<GNNEdgeSoA>;

#endif  // DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h
