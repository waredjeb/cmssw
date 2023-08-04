#ifndef IOMC_ParticleGun_CloseByParticleGunProducer_H
#define IOMC_ParticleGun_CloseByParticleGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

namespace edm {

  class CloseByParticleGunProducer : public BaseFlatGunProducer {
  public:
    CloseByParticleGunProducer(const ParameterSet&);
    ~CloseByParticleGunProducer() override;

    static void fillDescriptions(ConfigurationDescriptions& descriptions);

  private:
    void produce(Event& e, const EventSetup& es) override;

  protected:
    // data members
    bool fControlledByEta;
    double fEnMin, fEnMax, fEtaMin, fEtaMax, fRMin, fRMax, fZMin, fZMax, fDelta, fPhiMin, fPhiMax, fTMin, fTMax;
    int fNParticles;
    bool fMaxEnSpread = false;
    bool fPointing = false;
    bool fOverlapping = false;
    bool fRandomShoot = false;
    bool fUseDeltaT = true;
    double fOffsetFirst = 0.;
    std::vector<int> fPartIDs;
  };
}  // namespace edm

#endif
