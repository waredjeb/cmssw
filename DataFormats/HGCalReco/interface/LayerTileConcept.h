
#pragma once

namespace ticl::concepts {

  template <typename T>
  concept LayerTile = requires {
    T::nEtaBins;
    T::minEta;
    T::maxEta;
    T::nPhiBins;
    T::nBins;
  };

}
