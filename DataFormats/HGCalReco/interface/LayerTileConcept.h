
#pragma once

namespace ticl::concepts {

  template <typename T>
  concept LayerTile = requires {
    T::nEtaBins;
    T::minEta;
    T::maxEta;
    T::nPhiBins;
    T::nBins;
    T::nLayers;
    T::iterations;
  };

}  // namespace ticl::concepts
