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

      outTracksters.emplace_back(std::move(cleaned));
      outMap.emplace_back(std::move(kept));
      continue;
    }

    ticl::Trackster cleaned = link;
    std::vector<unsigned int> kept, dropped;
    std::vector<float> keptW; 
    kept.reserve(members.size());
    keptW.reserve(members.size());

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

    // compute energy of cleaned link
    double eNew = 0.0;
    if (!weightMode_) {
      for (auto idx : kept) eNew += in.clue3d[idx].raw_energy();
    } else {
      for (size_t i = 0; i < kept.size(); ++i)
        eNew += keptW[i] * in.clue3d[kept[i]].raw_energy();
    }
    setLinkRawEnergy_(cleaned, eNew);

    if (!weightMode_) keptW.clear();
    kept.shrink_to_fit();
    keptW.shrink_to_fit();

    outTracksters.emplace_back(std::move(cleaned));
    outMap.emplace_back(std::move(kept));

    if (emitDroppedAsStandalone_ && !dropped.empty()) {
      ticl::Trackster droppedLink = link;

      double eDrop = 0.0;
      for (auto idx : dropped) eDrop += in.clue3d[idx].raw_energy();
      setLinkRawEnergy_(droppedLink, eDrop);

      dropped.shrink_to_fit();

      outTracksters.emplace_back(std::move(droppedLink));
      outMap.emplace_back(std::move(dropped));
    }
  }
}
