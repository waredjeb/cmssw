#ifndef IOMC_ParticleGun_MultiCloseByParticleGunProducer_H
#define IOMC_ParticleGun_MultiCloseByParticleGunProducer_H

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

namespace edm {

  class MultiCloseByParticleGunProducer : public BaseFlatGunProducer {
  public:
    MultiCloseByParticleGunProducer(const ParameterSet&);
    ~MultiCloseByParticleGunProducer() override;

    static void fillDescriptions(ConfigurationDescriptions& descriptions);

  private:
    void produce(Event& e, const EventSetup& es) override;

  protected:
    // data members
    bool fControlledByEta;
    //define the following a std::vector<double> to allow for multiple particles fVarMin, fVarMax, fEtaMin, fEtaMax, fRMin, fRMax, fZMin, fZMax, fDelta, fPhiMin, fPhiMax, fTMin, fTMax, fOffsetFirst;
    std::vector<double> fVarMax;
    std::vector<double> fVarMin;
    std::vector<double> fEtaMax;
    std::vector<double> fEtaMin;
    std::vector<double> fRMax;
    std::vector<double> fRMin;
    std::vector<double> fZMax;
    std::vector<double> fZMin;
    double fDelta;
    std::vector<double> fPhiMin;
    std::vector<double> fPhiMax;
    std::vector<double> fTMin;
    std::vector<double> fTMax;
    std::vector<double> fOffsetFirst;
    int fNParticles;
    bool fMaxVarSpread = false;
    bool fFlatPtGeneration = false;
    bool fPointing = false;
    bool fOverlapping = false;
    bool fRandomShoot = false;
    bool fUseDeltaT = false;
    std::vector<int> fPartIDs;

    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_fieldToken;
  };
}  // namespace edm

#endif
