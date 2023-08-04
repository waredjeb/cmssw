#include <ostream>
#include <cmath>

#include "IOMC/ParticleGuns/interface/CloseByParticleGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

CloseByParticleGunProducer::CloseByParticleGunProducer(const ParameterSet& pset) : BaseFlatGunProducer(pset) {
  ParameterSet pgun_params = pset.getParameter<ParameterSet>("PGunParameters");
  fControlledByEta = pgun_params.getParameter<bool>("ControlledByEta");
  fEnMax = pgun_params.getParameter<double>("EnMax");
  fEnMin = pgun_params.getParameter<double>("EnMin");
  fMaxEnSpread = pgun_params.getParameter<bool>("MaxEnSpread");
  if (fControlledByEta) {
    fEtaMax = pgun_params.getParameter<double>("MaxEta");
    fEtaMin = pgun_params.getParameter<double>("MinEta");
    if (fEtaMax <= fEtaMin)
      LogError("CloseByParticleGunProducer") << " Please fix MinEta and MaxEta values in the configuration";
  } else {
    fRMax = pgun_params.getParameter<double>("RMax");
    fRMin = pgun_params.getParameter<double>("RMin");
    if (fRMax <= fRMin)
      LogError("CloseByParticleGunProducer") << " Please fix RMin and RMax values in the configuration";
  }
  fZMax = pgun_params.getParameter<double>("ZMax");
  fZMin = pgun_params.getParameter<double>("ZMin");
  fDelta = pgun_params.getParameter<double>("Delta");
  fPhiMin = pgun_params.getParameter<double>("MinPhi");
  fPhiMax = pgun_params.getParameter<double>("MaxPhi");
  fPointing = pgun_params.getParameter<bool>("Pointing");
  fOverlapping = pgun_params.getParameter<bool>("Overlapping");
  fRandomShoot = pgun_params.getParameter<bool>("RandomShoot");
  fNParticles = pgun_params.getParameter<int>("NParticles");
  fPartIDs = pgun_params.getParameter<vector<int> >("PartID");

  // set dt between particles
  fUseDeltaT = pgun_params.getParameter<bool>("UseDeltaT");
  fTMax = pgun_params.getParameter<double>("TMax");
  fTMin = pgun_params.getParameter<double>("TMin");
  fOffsetFirst = pgun_params.getParameter<double>("OffsetFirst");

  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

CloseByParticleGunProducer::~CloseByParticleGunProducer() {
  // no need to cleanup GenEvent memory - done in HepMCProduct
}

void CloseByParticleGunProducer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("AddAntiParticle", false);
  {
    edm::ParameterSetDescription psd0;
    psd0.add<bool>("ControlledByEta", false);
    psd0.add<double>("Delta", 10);
    psd0.add<double>("EnMax", 200.0);
    psd0.add<double>("EnMin", 25.0);
    psd0.add<bool>("MaxEnSpread", false);
    psd0.add<double>("MaxEta", 2.7);
    psd0.add<double>("MaxPhi", 3.14159265359);
    psd0.add<double>("MinEta", 1.7);
    psd0.add<double>("MinPhi", -3.14159265359);
    psd0.add<int>("NParticles", 2);
    psd0.add<bool>("Overlapping", false);
    psd0.add<std::vector<int>>("PartID", {
      22,
    });
    psd0.add<bool>("Pointing", true);
    psd0.add<double>("RMax", 120);
    psd0.add<double>("RMin", 60);
    psd0.add<bool>("RandomShoot", false);
    psd0.add<double>("ZMax", 321);
    psd0.add<double>("ZMin", 320);
    psd0.add<bool>("UseDeltaT", true);
    psd0.add<double>("TMin", 0.05);
    psd0.add<double>("TMax", 0.05);
    psd0.add<double>("OffsetFirst", 0.);
    desc.add<edm::ParameterSetDescription>("PGunParameters", psd0);
  }
  desc.addUntracked<int>("Verbosity", 0);
  desc.addUntracked<unsigned int>("firstRun", 1);
  desc.add<std::string>("psethack", "random particles in phi and r windows");
  descriptions.add("CloseByParticleGunProducer", desc);
  // or use the following to generate the label from the module's C++ type
  // descriptions.addWithDefaultLabel(desc);
}

void CloseByParticleGunProducer::produce(Event& e, const EventSetup& es) {
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  if (fVerbosity > 0) {
    LogDebug("CloseByParticleGunProducer") << " CloseByParticleGunProducer : Begin New Event Generation" << endl;
  }
  fEvt = new HepMC::GenEvent();

  int barcode = 1;
  unsigned int numParticles = fRandomShoot ? CLHEP::RandFlat::shoot(engine, 1, fNParticles) : fNParticles;

  double phi = CLHEP::RandFlat::shoot(engine, fPhiMin, fPhiMax);
  double fZ = CLHEP::RandFlat::shoot(engine, fZMin, fZMax);
  double fR;
  double fT;

  if (!fControlledByEta) {
    fR = CLHEP::RandFlat::shoot(engine, fRMin, fRMax);
  } else {
    double fEta = CLHEP::RandFlat::shoot(engine, fEtaMin, fEtaMax);
    fR = (fZ / sinh(fEta));
  }

  if (fUseDeltaT) {
    fT = CLHEP::RandFlat::shoot(engine, fTMin, fTMax);
  } else { 
    fT = 0.;
  }

  double tmpPhi = phi;
  double tmpR = fR;

  // Loop over particles
  for (unsigned int ip = 0; ip < numParticles; ++ip) {
    if (fOverlapping) {
      fR = CLHEP::RandFlat::shoot(engine, tmpR - fDelta, tmpR + fDelta);
      phi = CLHEP::RandFlat::shoot(engine, tmpPhi - fDelta / fR, tmpPhi + fDelta / fR);
    } else
      phi += fDelta / fR;

    double fEn;
    if (numParticles > 1 && fMaxEnSpread)
      fEn = fEnMin + ip * (fEnMax - fEnMin) / (numParticles - 1);
    else
      fEn = CLHEP::RandFlat::shoot(engine, fEnMin, fEnMax);

    int partIdx = CLHEP::RandFlat::shoot(engine, 0, fPartIDs.size());
    int PartID = fPartIDs[partIdx];
    const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID)));
    double mass = PData->mass().value();
    double mom2 = fEn * fEn - mass * mass;
    double mom = 0.;
    if (mom2 > 0.) {
      mom = sqrt(mom2);
    }
    double px = 0.;
    double py = 0.;
    double pz = mom;
    double energy = fEn;

    // Compute Vertex Position
    double x = fR * cos(phi);
    double y = fR * sin(phi);
    
    HepMC::FourVector p(px, py, pz, energy);
    // If we are requested to be pointing to (0,0,0), correct the momentum direction
    if (fPointing) {
      math::XYZVector direction(x, y, fZ);
      math::XYZVector momentum = direction.unit() * mom;
      p.setX(momentum.x());
      p.setY(momentum.y());
      p.setZ(momentum.z());
    }

    // compute correct path considering magnetic field 
    const double v = p.pz() / p.e() * c_light / cm;
    const double radius = sqrt(p.px() * p.px() + p.py() * p.py()) * 87.78f;  // cm (1 GeV track has 1 GeV/c / (e * 3.8T) ~ 87 cm radius in a 3.8T field)
    const double arc = 2 * asin(sqrt(x * x + y * y) / (2 * radius)) * radius;
    const double path = PData->charge() ? sqrt(arc * arc + fZ * fZ) : sqrt(x * x + y * y + fZ * fZ); 
    // if not pointing this doesn't mean a lot, keep the old way
    const double TimePath = fPointing ? (path / v) : (sqrt(x * x + y * y + fZ * fZ) / v); 
    double timeOffset = fOffsetFirst + (TimePath + ip * fT) * ns * c_light;
    // ns = 1, cm = 10, c_light is in mm/ns

    // std::cout << "shot particle from (" << x << ", " << y << ", " << fZ << "), with energy " << energy << ", at time " << TimePath << std::endl;
   
    HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(x * cm, y * cm, fZ * cm, timeOffset));

    HepMC::GenParticle* Part = new HepMC::GenParticle(p, PartID, 1);
    Part->suggest_barcode(barcode);
    barcode++;

    Vtx->add_particle_out(Part);

    if (fVerbosity > 0) {
      Vtx->print();
      Part->print();
    }
    fEvt->add_vertex(Vtx);
  }

  fEvt->set_event_number(e.id().event());
  fEvt->set_signal_process_id(20);

  if (fVerbosity > 0) {
    fEvt->print();
  }

  unique_ptr<HepMCProduct> BProduct(new HepMCProduct());
  BProduct->addHepMCData(fEvt);
  e.put(std::move(BProduct), "unsmeared");

  unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(std::move(genEventInfo));

  if (fVerbosity > 0) {
    LogDebug("CloseByParticleGunProducer") << " CloseByParticleGunProducer : Event Generation Done " << endl;
  }
}
