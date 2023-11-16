// Authors: Marco Rovere - marco.rovere@cern.ch, Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 11/2018

#ifndef RecoHGCal_TICL_ClusterFilterByAlgoAndSize_H__
#define RecoHGCal_TICL_ClusterFilterByAlgoAndSize_H__

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "PreliminaryClusterFilterBase.h"

#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class PreliminaryClusterFilterByCP : public PreliminaryClusterFilterBase {
  public:
    PreliminaryClusterFilterByCP (const edm::ParameterSet& ps)
        : PreliminaryClusterFilterBase(ps),
        delta_eta_(ps.getParameter<double>("deltaEta")),
        delta_phi_(ps.getParameter<double>("deltaPhi")){}

    ~PreliminaryClusterFilterByCP () override{};

    void filter(const std::vector<reco::CaloCluster>& layerClusters,
                std::vector<float>& layerClustersMask,
                const std::vector<CaloParticle>& caloparticles,
                const std::vector<size_t>& caloParticlesIndices) const override {
      std::array<TICLLayerTile, 2> lcTiles = {};      // all Tracksters, propagated to lastLayerEE
      for (size_t iCL = 0; iCL < layerClustersMask.size(); ++iCL) {
        auto const& cl = layerClusters[iCL];
        lcTiles[(cl.z() > 0 ? 1 : 0)].fill(cl.eta(), cl.phi(), iCL);   
      }
      for (auto const& cpIndex : caloParticlesIndices) {
        auto const& cp = caloparticles[cpIndex];
        auto const etaCP = cp.eta(); 
        auto const phiCP = cp.phi(); 
        auto const& tile = lcTiles[etaCP > 0]; 
        float eta_min = std::max(abs(etaCP) - delta_eta_, (float)TileConstants::minEta);
        float eta_max = std::min(abs(etaCP) + delta_eta_, (float)TileConstants::maxEta);
        std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, eta_max, phiCP - delta_phi_, phiCP + delta_phi_);
        for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
          for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
            const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
            for (const unsigned &t_i : in_tile) {
              auto const& cl = layerClusters[t_i];
              layerClustersMask[t_i] = 1.f;
            }
          }
        }
      }
    }

  private:
    float delta_eta_;
    float delta_phi_;
  };
}  // namespace ticl

#endif
