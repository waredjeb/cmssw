#ifndef DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h
#define DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/TrackstersSoA.h"

using TrackstersSoAHostCollection = PortableMultiCollection<alpaka::DevCpu, GNNNodeSoA, GNNEdgeSoA, GNNEdgeIndexSoA>;
using TrackstersGNNOutputSoAHostCollection = PortableHostCollection<GNNOutputSoA>;
using TrackstersGNNPostprocessingSoAHostCollection = PortableHostCollection<GNNPostprocessingSoA>;

#endif  // DataFormats_HGCalReco_interface_TrackstersSoAHostCollection_h
