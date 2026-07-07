#ifndef DataFormats_PortableTestObjects_interface_alpaka_TrackstersSoADeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_TrackstersSoADeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/TrackstersSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using TrackstersSoADeviceCollection = PortableCollection<GNNTrackstersBlocks>;
  using TrackstersGNNOutputSoADeviceCollection = PortableCollection<GNNOutputSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_PortableTestObjects_interface_alpaka_TrackstersSoADeviceCollection_h
