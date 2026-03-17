#ifndef DataFormats_HGCalReco_interface_alpaka_CLUE3DStateDeviceCollection_h
#define DataFormats_HGCalReco_interface_alpaka_CLUE3DStateDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/HGCalReco/interface/CLUE3DStateSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using CLUE3DStateDeviceCollection = PortableCollection<CLUE3DStateSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_HGCalReco_interface_alpaka_CLUE3DStateDeviceCollection_h
