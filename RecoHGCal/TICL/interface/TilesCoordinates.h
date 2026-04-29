
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <vector>

namespace ticl {

  struct TilesCoordinates {
    std::vector<float> etas_pos;
    std::vector<float> etas_neg;
    std::vector<float> phis_pos;
    std::vector<float> phis_neg;
    std::vector<std::uint32_t> ids_pos;
    std::vector<std::uint32_t> ids_neg;

    TilesCoordinates(std::integral auto size) {
      etas_pos.reserve(size);
      etas_neg.reserve(size);
      phis_pos.reserve(size);
      phis_neg.reserve(size);
      ids_pos.reserve(size);
      ids_neg.reserve(size);
    }

    auto size() const { return std::array<std::size_t, 2>{etas_neg.size(), etas_pos.size()}; }
  };

}  // namespace ticl
