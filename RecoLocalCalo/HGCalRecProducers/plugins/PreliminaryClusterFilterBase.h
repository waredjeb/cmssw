// Author: Felice Pantaleo, Marco Rovere - felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 09/2018

#ifndef RecoHGCal_TICL_PreliminaryClusterFilterBase_H__
#define RecoHGCal_TICL_PreliminaryClusterFilterBase_H__

#include "DataFormats/HGCalReco/interface/Common.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include <memory>
#include <vector>

namespace edm {
  class ParameterSet;
}
namespace reco {
  class CaloCluster;
}

namespace ticl {
  class PreliminaryClusterFilterBase {
  public:
    explicit PreliminaryClusterFilterBase(const edm::ParameterSet&){};
    virtual ~PreliminaryClusterFilterBase(){};

    virtual void filter(const std::vector<reco::CaloCluster>& layerClusters,
                        std::vector<float>& layerClustersMask,
                        const std::vector<CaloParticle>& caloparticles,
                        const std::vector<size_t>& caloParticlesIndices) const = 0;
  };
}  // namespace ticl

#endif
