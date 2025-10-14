// Author: Felice Pantaleo (CERN), 2025, felice.pantaleo@cern.ch
#ifndef DataFormats_Common_MultiCollection_h
#define DataFormats_Common_MultiCollection_h

#include <span>
#include <utility>
#include <vector>

#include "DataFormats/Common/interface/RefProd.h"
#include "MultiSpan.h"

namespace edm {

/**
 * @brief Lightweight persistent holder for several `edm::RefProd<Collection>`
 *        objects.
 *
 * Only the vector of `RefProd`s is stored on disk.  No transient caches or
 * synchronisation primitives live inside the class, so it is trivially
 * movable.  Consumers obtain a flat view with `makeFlatView()`, which assembles
 * and *returns* a fully‑populated `MultiSpan` by value.
 *
 * This design avoids copy/move issues with `std::once_flag` and makes the type
 * EDM‑wrapper‑friendly while still giving fast, indexed access to the
 * concatenated elements.
 */
template <typename Collection>
class MultiCollection {
public:
  using value_type = typename Collection::value_type;

  MultiCollection() = default;

  explicit MultiCollection(std::initializer_list<edm::RefProd<Collection>> refs) : refProds_{refs} {}

  // ---------------- producer‑side API ----------------------------------
  void add(edm::RefProd<Collection> const& ref) { refProds_.push_back(ref); }

  // ---------------- consumer‑side helpers ------------------------------
  /**
   * @brief Build and return a flat view that spans all referenced collections.
   *
   * The returned `MultiSpan` is independent of `this`, so callers may
   * move or store it locally without keeping the manager alive.
   */
  [[nodiscard]] MultiSpan<value_type> makeFlatView() const {
    MultiSpan<value_type> ms;
    for (auto const& rp : refProds_) {
      auto const& coll = *rp;  // Framework‑managed retrieval
      ms.add(std::span<const value_type>(coll.data(), coll.size()));
    }
    return ms;
  }

  const std::vector<edm::RefProd<Collection>>& refProds() const { return refProds_; }

private:
  std::vector<edm::RefProd<Collection>> refProds_;
};

}  // namespace edm

#endif  // DataFormats_Common_MultiCollection_h
