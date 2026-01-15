#include <vector>
#include <cmath>
#include <algorithm>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoHGCal/TICL/plugins/TracksterCleaningByBeta.h"
#include "RecoHGCal/TICL/interface/TracksterCleaningAlgoBase.h"

namespace {
  constexpr double c_cm_per_ns     = 29.9792458;
  constexpr double inv_c_cm_per_ns = 1.0 / c_cm_per_ns;
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
      drPower_prune_(conf.getParameter<double>("drPower_prune"))
{}

void TracksterCleaningByBeta::cleanTracksters(const Inputs& in,
                                              std::vector<ticl::Trackster>& outTracksters,
                                              std::vector<std::vector<unsigned int>>& outMap) const {
  const size_t nL = in.linked.size();
  outTracksters.clear(); outMap.clear();
  outTracksters.reserve(nL * (emitDroppedAsStandalone_ ? 2 : 1));
  outMap.reserve(nL * (emitDroppedAsStandalone_ ? 2 : 1));

  for (size_t L = 0; L < nL; ++L) {
    const auto& link    = in.linked[L];
    const auto& members = in.map[L];

    const auto& bcL = link.barycenter();
    const double LL = std::sqrt(bcL.x()*bcL.x() + bcL.y()*bcL.y() + bcL.z()*bcL.z());
    const double tLcorr = link.time() - LL * inv_c_cm_per_ns;
    const double zL  = bcL.z();
    const double etaL = bcL.eta();
    const double phiL = bcL.phi();

    // beta = log( E_link * sum(deltaR) )
    const double Ek = std::max(epsE_, linkEnergy_(link));
    double sumDR = 0.0;
    const size_t nm = members.size();
    if (nm >= 2) {
      for (size_t a = 0; a + 1 < nm; ++a) {
        const auto& ta = in.clue3d[members[a]];
        const double etaA = ta.barycenter().eta();
        const double phiA = ta.barycenter().phi();
        for (size_t b = a + 1; b < nm; ++b) {
          const auto& tb = in.clue3d[members[b]];
          const double dR = reco::deltaR(etaA, phiA, tb.barycenter().eta(), tb.barycenter().phi());
          if (dR <= R0_) sumDR += dR;
        }
      }
    }
    const double beta = std::log(Ek * std::max(epsDR_, sumDR));

    if (beta < betaContamMin_) {
      std::vector<unsigned int> kept = members;
      kept.shrink_to_fit();
      ticl::Trackster cleaned = link;

      // cleaned vertices to match CLUE3D membership
      cleaned.vertices().assign(members.begin(), members.end());
      cleaned.vertex_multiplicity().assign(members.size(), 1.0f);

      // clear edges
      cleaned.edges().clear();

      cleaned.calculateRawPt();
      cleaned.calculateRawEmPt();

      outTracksters.emplace_back(std::move(cleaned));
      outMap.emplace_back(std::move(kept));
      continue;
    }

    std::vector<unsigned int> kept, dropped;
    std::vector<float> keptW;
    kept.reserve(members.size());
    keptW.reserve(members.size());
    dropped.reserve(members.size());

    const double sZ_w = std::max(1e-12, sigmaZ_);
    const double sT_w = std::max(1e-12, sigmaT_);
    const double sR_w = std::max(1e-12, sigmaDR_);
    const double zPw = zPower_, tPw = tPower_, rPw = drPower_;

    const double sZ_p = std::max(1e-12, pruneUseSeparateKernels_ ? sigmaZ_prune_ : sigmaZ_);
    const double sT_p = std::max(1e-12, pruneUseSeparateKernels_ ? sigmaT_prune_ : sigmaT_);
    const double sR_p = std::max(1e-12, pruneUseSeparateKernels_ ? sigmaDR_prune_ : sigmaDR_);
    const double zPp = pruneUseSeparateKernels_ ? zPower_prune_ : zPower_;
    const double tPp = pruneUseSeparateKernels_ ? tPower_prune_ : tPower_;
    const double rPp = pruneUseSeparateKernels_ ? drPower_prune_ : drPower_;

    for (const auto idx : members) {
      const auto& t  = in.clue3d[idx];
      const auto& bc = t.barycenter();

      const double dz = bc.z() - zL;

      const double Lm = std::sqrt(bc.x()*bc.x() + bc.y()*bc.y() + bc.z()*bc.z());
      const double tmCorr = t.time() - Lm * inv_c_cm_per_ns;
      const double dt = tmCorr - tLcorr;

      const double dR = reco::deltaR(etaL, phiL, bc.eta(), bc.phi());

      const bool passZ = std::abs(dz) <= zAbsCut_;
      const bool passT = (tAbsCut_ <= 0.0) ? true : (std::abs(dt) <= tAbsCut_);

      // pruning 
      if (doPruning_) {
        double wzp = std::exp(-0.5 * (dz*dz)/(sZ_p*sZ_p));
        double wtp = std::exp(-0.5 * (dt*dt)/(sT_p*sT_p));
        double wrp = std::exp(-0.5 * (dR*dR)/(sR_p*sR_p));
        if (!passZ) wzp = 0.0;
        if (!passT) wtp = 0.0;
        const double w_prune = std::pow(wzp, zPp) * std::pow(wtp, tPp) * std::pow(wrp, rPp);
        if (w_prune < pruneWmin_) {
          dropped.push_back(idx);
          continue;
        }
      }
      
      // weighting
      double wzw = std::exp(-0.5 * (dz*dz)/(sZ_w*sZ_w));
      double wtw = std::exp(-0.5 * (dt*dt)/(sT_w*sT_w));
      double wrw = std::exp(-0.5 * (dR*dR)/(sR_w*sR_w));
      if (!passZ) wzw = 0.0;
      if (!passT) wtw = 0.0;
      const double w_energy = std::pow(wzw, zPw) * std::pow(wtw, tPw) * std::pow(wrw, rPw);

      if (!doPruning_) {
        if (w_energy >= wmin_) {
          kept.push_back(idx);
          if (weightMode_) keptW.push_back(static_cast<float>(w_energy));
        } else {
          dropped.push_back(idx);
        }
      } else {
        kept.push_back(idx);
        if (weightMode_) keptW.push_back(static_cast<float>(w_energy));
      }
    }

    if (kept.empty()) {
      kept = members;
      keptW.clear();
    }

    ticl::Trackster cleaned = link;

    cleaned.vertices().assign(kept.begin(), kept.end());
    cleaned.vertex_multiplicity().assign(kept.size(), 1.0f);

    cleaned.edges().clear();

    // compute energy of cleaned link
    double eNew = 0.0;

    double wx = 0.0, wy = 0.0, wz = 0.0;
    double wsum_pos = 0.0;

    double wt = 0.0;
    double wsum_t = 0.0;

    auto weight_i = [&](size_t i)->double {
      if (!weightMode_) return 1.0;
      return (i < keptW.size() ? static_cast<double>(keptW[i]) : 1.0);
    };

    for (size_t i = 0; i < kept.size(); ++i) {
      const unsigned int idx = kept[i];
      const auto& ts = in.clue3d[idx];
      const auto& bc = ts.barycenter();

      const double w = weight_i(i);
      const double e = static_cast<double>(ts.raw_energy());

      // reweighted energy
      const double e_eff = weightMode_ ? (w * e) : e;
      eNew += e_eff;

      // barycenter
      wx += e_eff * bc.x();
      wy += e_eff * bc.y();
      wz += e_eff * bc.z();
      wsum_pos += e_eff;

      // time
      wt += e_eff * ts.time();
      wsum_t += e_eff;
    }

    if (wsum_pos > 0.0) {
      cleaned.setBarycenter(ticl::Trackster::Vector(wx/wsum_pos, wy/wsum_pos, wz/wsum_pos));
    } else {
      cleaned.setBarycenter(link.barycenter());
    }

    if (wsum_t > 0.0) {
      cleaned.setTimeAndError(static_cast<float>(wt/wsum_t), -1.f);
    }

    setLinkRawEnergy_(cleaned, eNew);
    cleaned.calculateRawPt();
    cleaned.calculateRawEmPt();
    cleaned.zeroProbabilities(); 

    // store outputs
    outTracksters.emplace_back(std::move(cleaned));
    outMap.emplace_back(std::move(kept));

    // optionally emit dropped as a standalone link
    if (emitDroppedAsStandalone_ && !dropped.empty()) {
      ticl::Trackster droppedLink = link;

      droppedLink.vertices().assign(dropped.begin(), dropped.end());
      droppedLink.vertex_multiplicity().assign(dropped.size(), 1.0f);
      droppedLink.edges().clear();
      droppedLink.zeroProbabilities();

      double eDrop = 0.0;
      double dx = 0.0, dy = 0.0, dz = 0.0, dsum = 0.0;
      double dt = 0.0, dtsum = 0.0;
      for (auto idx : dropped) {
        const auto& ts = in.clue3d[idx];
        const auto& bc = ts.barycenter();
        const double e = static_cast<double>(ts.raw_energy());
        eDrop += e;
        dx += e * bc.x();
        dy += e * bc.y();
        dz += e * bc.z();
        dsum += e;
        dt += e * ts.time();
        dtsum += e;
      }
      if (dsum > 0.0) {
        droppedLink.setBarycenter(ticl::Trackster::Vector(dx/dsum, dy/dsum, dz/dsum));
      }
      if (dtsum > 0.0) {
        droppedLink.setTimeAndError(static_cast<float>(dt/dtsum), -1.f);
      }

      setLinkRawEnergy_(droppedLink, eDrop);
      droppedLink.calculateRawPt();
      droppedLink.calculateRawEmPt();

      outTracksters.emplace_back(std::move(droppedLink));
      outMap.emplace_back(std::move(dropped));
    }
  }
}
