#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

#include "DataFormats/Math/interface/deltaR.h"

#include "RecoHGCal/TICL/plugins/TracksterCleaningByBeta.h"
#include "RecoHGCal/TICL/interface/TracksterCleaningAlgoBase.h"

namespace {
  constexpr double c_cm_per_ns     = 29.9792458;
  constexpr double inv_c_cm_per_ns = 1.0 / c_cm_per_ns;

  inline bool validTime(double t) {
    return std::isfinite(t) && t > -90.0;
  }
}

using namespace ticl;

TracksterCleaningByBeta::TracksterCleaningByBeta(const edm::ParameterSet& conf, edm::ConsumesCollector iC)
    : TracksterCleaningAlgoBase(conf, iC),
      betaContamMin_(conf.getParameter<double>("betaContamMin")),
      R0_(conf.getParameter<double>("R0")),
      epsE_(conf.getParameter<double>("epsE")),
      epsDR_(conf.getParameter<double>("epsDR")),
      useRawEnergy_(conf.getParameter<bool>("useRawEnergy")),
      emitDroppedAsStandalone_(conf.getParameter<bool>("emitDroppedAsStandalone")),
      weightMode_(conf.getParameter<bool>("weightMode")),
      zAbsCut_(conf.getParameter<double>("zAbsCut")),
      tAbsCut_(conf.getParameter<double>("tAbsCut")),
      sigmaZ_(conf.getParameter<double>("sigmaZ")),
      sigmaT_(conf.getParameter<double>("sigmaT")),
      sigmaDR_(conf.getParameter<double>("sigmaDR")),
      zPower_(conf.getParameter<double>("zPower")),
      tPower_(conf.getParameter<double>("tPower")),
      drPower_(conf.getParameter<double>("drPower")),
      wmin_(conf.getParameter<double>("wmin")),
      doPruning_(conf.getParameter<bool>("doPruning")),
      pruneWmin_(conf.getParameter<double>("pruneWmin")),
      pruneUseSeparateKernels_(conf.getParameter<bool>("pruneUseSeparateKernels")),
      sigmaZ_prune_(conf.getParameter<double>("sigmaZ_prune")),
      sigmaT_prune_(conf.getParameter<double>("sigmaT_prune")),
      sigmaDR_prune_(conf.getParameter<double>("sigmaDR_prune")),
      zPower_prune_(conf.getParameter<double>("zPower_prune")),
      tPower_prune_(conf.getParameter<double>("tPower_prune")),
      drPower_prune_(conf.getParameter<double>("drPower_prune")) {}

void TracksterCleaningByBeta::cleanTracksters(const Inputs& in,
                                              std::vector<ticl::Trackster>& outTracksters,
                                              std::vector<std::vector<unsigned int>>& outMap) const {
  auto rebuildFromMembers =
      [&](ticl::Trackster const& seed,
          std::vector<unsigned int> const& keptMembers,
          std::vector<float> const* keptMemberWeights) -> ticl::Trackster {
        if (keptMembers.empty()) return seed;

        ticl::Trackster out = seed;
        out.vertices().clear();
        out.vertex_multiplicity().clear();
        out.edges().clear();

        const auto nLC = in.layerClusters.size();

        thread_local std::vector<float> lcWeight;
        thread_local std::vector<unsigned int> touchedLCs;

        if (lcWeight.size() != nLC) {
          lcWeight.assign(nLC, 0.f);
          touchedLCs.clear();
        } else {
          for (auto idx : touchedLCs) lcWeight[idx] = 0.f;
          touchedLCs.clear();
        }

        for (size_t im = 0; im < keptMembers.size(); ++im) {
          const unsigned int mIdx = keptMembers[im];
          if (mIdx >= in.clue3d.size()) continue;

          auto const& ts = in.clue3d[mIdx];
          const float wMember =
              (keptMemberWeights && im < keptMemberWeights->size()) ? (*keptMemberWeights)[im] : 1.f;

          for (auto lcIdx : ts.vertices()) {
            if (lcIdx >= nLC) continue;
            if (lcWeight[lcIdx] == 0.f) touchedLCs.push_back(lcIdx);
            if (wMember > lcWeight[lcIdx]) lcWeight[lcIdx] = wMember;
          }
        }

        if (touchedLCs.empty()) return seed;

        std::sort(touchedLCs.begin(), touchedLCs.end());
        touchedLCs.erase(std::unique(touchedLCs.begin(), touchedLCs.end()), touchedLCs.end());

        out.vertices().reserve(touchedLCs.size());
        out.vertex_multiplicity().reserve(touchedLCs.size());

        double eSum = 0.0;
        double wx = 0.0, wy = 0.0, wz = 0.0;

        for (auto lcIdx : touchedLCs) {
          auto const& lc  = in.layerClusters[lcIdx];
          auto const& pos = lc.position();

          const double w   = double(lcWeight[lcIdx]);
          const double eLC = double(lc.energy()) * w;

          out.vertices().push_back(lcIdx);
          out.vertex_multiplicity().push_back(1.f);

          eSum += eLC;
          wx += eLC * pos.x();
          wy += eLC * pos.y();
          wz += eLC * pos.z();
        }

        if (eSum > 0.0)
          out.setBarycenter(ticl::Trackster::Vector(wx / eSum, wy / eSum, wz / eSum));
        else
          out.setBarycenter(seed.barycenter());

        setLinkRawEnergy_(out, eSum);
        out.calculateRawPt();
        out.zeroProbabilities();
        return out;
      };

  // You said: we're not doing pruning. Enforce it here so the behavior is unambiguous.
  if (doPruning_) {
    // If you ever accidentally turn it on in config, fail loudly.
    throw cms::Exception("TracksterCleaningByBeta")
        << "doPruning_=True but this debug build assumes pruning is disabled.";
  }

  const size_t nL = in.linked.size();
  outTracksters.clear();
  outMap.clear();
  outTracksters.reserve(nL * (emitDroppedAsStandalone_ ? 2 : 1));
  outMap.reserve(nL * (emitDroppedAsStandalone_ ? 2 : 1));

  for (size_t L = 0; L < nL; ++L) {
    auto const& link    = in.linked[L];
    auto const& members = in.map[L];

    // Stay safe.
    if (members.empty()) {
      outTracksters.emplace_back(link);
      outMap.emplace_back(members);
      continue;
    }

    // Precompute link kinematics/time correction
    const auto& bcL = link.barycenter();
    const double LL = std::sqrt(bcL.x()*bcL.x() + bcL.y()*bcL.y() + bcL.z()*bcL.z());
    const double etaL = bcL.eta();
    const double phiL = bcL.phi();
    const double zL   = bcL.z();

    const bool linkHasT = validTime(link.time());
    const double tLcorr = linkHasT ? (link.time() - LL * inv_c_cm_per_ns) : 0.0;

    // Sum member energy (informational)
    double eMembersIn = 0.0;
    for (auto idx : members) {
      if (idx >= in.clue3d.size()) continue;
      eMembersIn += in.clue3d[idx].raw_energy();
    }

    // beta = log( E_link * sum(deltaR) )
    const double Ek = std::max(epsE_, linkEnergy_(link));
    double sumDR = 0.0;
    const size_t nm = members.size();
    if (nm >= 2) {
      for (size_t a = 0; a + 1 < nm; ++a) {
        const auto ia = members[a];
        if (ia >= in.clue3d.size()) continue;
        const auto& ta = in.clue3d[ia];
        const double etaA = ta.barycenter().eta();
        const double phiA = ta.barycenter().phi();
        for (size_t b = a + 1; b < nm; ++b) {
          const auto ib = members[b];
          if (ib >= in.clue3d.size()) continue;
          const auto& tb = in.clue3d[ib];
          const double dR = reco::deltaR(etaA, phiA, tb.barycenter().eta(), tb.barycenter().phi());
          if (dR <= R0_) sumDR += dR;
        }
      }
    }
    const double beta = std::log(Ek * std::max(epsDR_, sumDR));

    // beta >= betaMin => apply weight-based thresholding (no pruning)
    std::vector<unsigned int> kept, dropped;
    std::vector<float> keptW;

    kept.reserve(members.size());
    dropped.reserve(members.size());
    if (weightMode_) keptW.reserve(members.size());

    const double sZ = std::max(1e-12, sigmaZ_);
    const double sT = std::max(1e-12, sigmaT_);
    const double sR = std::max(1e-12, sigmaDR_);
    const double zP = zPower_, tP = tPower_, rP = drPower_;

    for (auto idx : members) {
      if (idx >= in.clue3d.size()) continue;
      const auto& t  = in.clue3d[idx];
      const auto& bc = t.barycenter();

      const double dz = bc.z() - zL;
      const double dR = reco::deltaR(etaL, phiL, bc.eta(), bc.phi());

      const bool passZ = (std::abs(dz) <= zAbsCut_);

      const bool memberHasT = validTime(t.time());
      const bool useTime = (tAbsCut_ > 0.0) && linkHasT && memberHasT;

      bool passT = true;
      double dt = 0.0;
      double wtw = 1.0;

      if (useTime) {
        const double Lm = std::sqrt(bc.x()*bc.x() + bc.y()*bc.y() + bc.z()*bc.z());
        const double tmCorr = t.time() - Lm * inv_c_cm_per_ns;
        dt = tmCorr - tLcorr;
        passT = (std::abs(dt) <= tAbsCut_);
        wtw = std::exp(-0.5 * (dt * dt) / (sT * sT));
      }

      double wzw = std::exp(-0.5 * (dz * dz) / (sZ * sZ));
      double wrw = std::exp(-0.5 * (dR * dR) / (sR * sR));
      if (!passZ) wzw = 0.0;
      if (!passT) wtw = 0.0;

      const double w_energy = std::pow(wzw, zP) * std::pow(wtw, tP) * std::pow(wrw, rP);

      const bool keep = (w_energy >= wmin_);

      if (keep) {
        kept.push_back(idx);
        if (weightMode_) keptW.push_back(static_cast<float>(w_energy));
      } else {
        dropped.push_back(idx);
      }
    }

    // Safety: if everything dropped, revert to keeping all (as original)
    if (kept.empty()) {
      kept = members;
      keptW.clear();
    }

    ticl::Trackster cleaned = rebuildFromMembers(link, kept, weightMode_ ? &keptW : nullptr);

    const double eBefore = link.raw_energy();
    const double eAfter  = cleaned.raw_energy();

    outTracksters.emplace_back(std::move(cleaned));
    outMap.emplace_back(std::move(kept));
    

    if (emitDroppedAsStandalone_ && !dropped.empty()) {
      ticl::Trackster droppedLink = rebuildFromMembers(link, dropped, nullptr);
      outTracksters.emplace_back(std::move(droppedLink));
      outMap.emplace_back(std::move(dropped));
    }
  }
}
