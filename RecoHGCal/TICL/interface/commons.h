#ifndef RecoHGCal_TICL_interface_commons_h
#define RecoHGCal_TICL_interface_commons_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {

  //constants
  constexpr double mpion = 0.13957;
  constexpr float mpion2 = mpion * mpion;
  enum LayerType {

    CE_E_120 = 0,
    CE_E_200 = 1,
    CE_E_300 = 2,
    CE_H_120_F = 3,
    CE_H_200_F = 4,
    CE_H_300_F = 5,
    CE_H_120_C = 6,
    CE_H_200_C = 7,
    CE_H_SCINT_C = 8,
    EnumSize = 9

  };

  inline Trackster::ParticleType tracksterParticleTypeFromPdgId(int pdgId, int charge) {
    if (pdgId == 111) {
      return Trackster::ParticleType::neutral_pion;
    } else {
      pdgId = std::abs(pdgId);
      if (pdgId == 22) {
        return Trackster::ParticleType::photon;
      } else if (pdgId == 11) {
        return Trackster::ParticleType::electron;
      } else if (pdgId == 13) {
        return Trackster::ParticleType::muon;
      } else {
        bool isHadron = (pdgId > 100 and pdgId < 900) or (pdgId > 1000 and pdgId < 9000);
        if (isHadron) {
          if (charge != 0) {
            return Trackster::ParticleType::charged_hadron;
          } else {
            return Trackster::ParticleType::neutral_hadron;
          }
        } else {
          return Trackster::ParticleType::unknown;
        }
      }
    }
  }

}  // namespace ticl

#endif
