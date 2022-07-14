#if !defined(SimCalorimetry_HGCalAssociatorProducers_interface_AssociatorTools_h)
#define SimCalorimetry_HGCalAssociatorProducers_interface_AssociatorTools_h
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include <vector>

static void removeCPFromPU(const std::vector<CaloParticle>& caloParticles,
                           std::vector<size_t>& cPIndices,
                           bool hardScatterOnly = true) {
  //Consider CaloParticles coming from the hard scatterer
  //excluding the PU contribution and save the indices.
  for (unsigned int cpId = 0; cpId < caloParticles.size(); ++cpId) {
    if (hardScatterOnly && (caloParticles[cpId].g4Tracks()[0].eventId().event() != 0 or
                            caloParticles[cpId].g4Tracks()[0].eventId().bunchCrossing() != 0)) {
      LogDebug("HGCalValidator") << "Excluding CaloParticles from event: "
                                 << caloParticles[cpId].g4Tracks()[0].eventId().event()
                                 << " with BX: " << caloParticles[cpId].g4Tracks()[0].eventId().bunchCrossing()
                                 << std::endl;
      continue;
    }
    cPIndices.emplace_back(cpId);
  }
}

namespace hgcal {

  // This introduces a CaloParticle on layer concept. For a CaloParticle it stores:
  // 1. Its id: caloParticleId.
  // 2. The energy that the CaloParticle deposited in a specific layer and it was reconstructed.
  // 3. The hits_and_fractions that contributed to that deposition. SimHits that aren't reconstructed
  //    and doesn't have any matched rechits are disregarded. Keep in mind that since a CaloParticle
  //    should most probably have more than one SimCluster, all different contributions from the same CaloParticle
  //    to a single hit are merged into a single entry, with the fractions properly summed.
  // 4. A map to save the LayerClusters ids (id is the key) that reconstructed at least one SimHit of the CaloParticle under study
  //    together with the energy that the LayerCluster reconstructed from the CaloParticle and the score. The energy
  //    is not the energy of the LayerCluster, but the energy of the LayerCluster coming from the CaloParticle.
  //    So, there will be energy of the LayerCluster that is disregarded here, since there may be LayerCluster's
  //    cells that the CaloParticle didn't contribute.
  struct simObjectOnLayer {
    unsigned int simObjectId;
    float energy = 0;
    std::vector<std::pair<DetId, float>> hits_and_fractions;
    std::unordered_map<int, std::pair<float, float>> clusterIdToEnergyAndScore;
  };

  typedef struct simObjectOnLayer caloParticleOnLayer;
  typedef struct simObjectOnLayer simClusterOnLayer;

}  // namespace hgcal
#endif
