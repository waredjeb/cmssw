#include <ostream>
#include <cmath>

#include "IOMC/ParticleGuns/interface/MultiCloseByParticleGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/GlobalPhysicalConstants.h>
#include <CLHEP/Random/RandFlat.h>

using namespace edm;
using namespace std;

MultiCloseByParticleGunProducer::MultiCloseByParticleGunProducer(const ParameterSet& pset)
    : BaseFlatGunProducer(pset), m_fieldToken(esConsumes()) {
  ParameterSet pgun_params = pset.getParameter<ParameterSet>("PGunParameters");
  fControlledByEta = pgun_params.getParameter<bool>("ControlledByEta");
  fVarMax = pgun_params.getParameter<std::vector<double>>("VarMax");
  fVarMin = pgun_params.getParameter<std::vector<double>>("VarMin");
  fMaxVarSpread = pgun_params.getParameter<bool>("MaxVarSpread");
  fFlatPtGeneration = pgun_params.getParameter<bool>("FlatPtGeneration");

  if (fControlledByEta) {
    fEtaMax = pgun_params.getParameter<std::vector<double>>("MaxEta");
    fEtaMin = pgun_params.getParameter<std::vector<double>>("MinEta");
    if (fEtaMax <= fEtaMin)
      LogError("MultiCloseByParticleGunProducer") << " Please fix MinEta and MaxEta values in the configuration";
  } else {
    fRMax = pgun_params.getParameter<std::vector<double>>("RMax");
    fRMin = pgun_params.getParameter<std::vector<double>>("RMin");
    if (fRMax <= fRMin)
      LogError("MultiCloseByParticleGunProducer") << " Please fix RMin and RMax values in the configuration";
  }
  fZMax = pgun_params.getParameter<std::vector<double>>("ZMax");
  fZMin = pgun_params.getParameter<std::vector<double>>("ZMin");
  fDelta = pgun_params.getParameter<double>("Delta");
  fPhiMin = pgun_params.getParameter<std::vector<double>>("MinPhi");
  fPhiMax = pgun_params.getParameter<std::vector<double>>("MaxPhi");
  fPointing = pgun_params.getParameter<bool>("Pointing");
  fOverlapping = pgun_params.getParameter<bool>("Overlapping");
  if (fFlatPtGeneration && !fPointing)
    LogError("MultiCloseByParticleGunProducer")
        << " Can't generate non pointing FlatPt samples; please disable FlatPt generation or generate pointing sample";
  fRandomShoot = pgun_params.getParameter<bool>("RandomShoot");
  fNParticles = pgun_params.getParameter<int>("NParticles");
  fPartIDs = pgun_params.getParameter<vector<int>>("PartID");
  for(size_t ip = 0; ip < fPartIDs.size(); ip++)
    if (fVarMin[ip] < 1 && !fFlatPtGeneration){
        LogError("MultiCloseByParticleGunProducer") << " Please choose a minimum energy greater than 1 GeV, otherwise time "
                                                      "information may be invalid or not reliable";
    }
  // set dt between particles
  fUseDeltaT = pgun_params.getParameter<bool>("UseDeltaT");
  fTMax = pgun_params.getParameter<std::vector<double>>("TMax");
  fTMin = pgun_params.getParameter<std::vector<double>>("TMin");
  if (fTMax <= fTMin)
    LogError("MultiCloseByParticleGunProducer") << " Please fix TMin and TMax values in the configuration";
  // set a fixed time offset for the particles
  fOffsetFirst = pgun_params.getParameter<std::vector<double>>("OffsetFirst");

  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

MultiCloseByParticleGunProducer::~MultiCloseByParticleGunProducer() {
  // no need to cleanup GenEvent memory - done in HepMCProduct
}

void MultiCloseByParticleGunProducer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("AddAntiParticle", false);
  {
    edm::ParameterSetDescription psd0;
    psd0.add<bool>("ControlledByEta", false);
    psd0.add<double>("Delta", 10);
    psd0.add<std::vector<double>>("VarMax", {200.0});
    psd0.add<std::vector<double>>("VarMin", {25.0});
    psd0.add<bool>("MaxVarSpread", false);
    psd0.add<bool>("FlatPtGeneration", false);
    psd0.add<std::vector<double>>("MaxEta", {2.7});
    psd0.add<std::vector<double>>("MaxPhi", {3.14159265359});
    psd0.add<std::vector<double>>("MinEta", {1.7});
    psd0.add<std::vector<double>>("MinPhi", {-3.14159265359});
    psd0.add<int>("NParticles", 2);
    psd0.add<bool>("Overlapping", false);
    psd0.add<std::vector<int>>("PartID",
                               {
                                   22,
                               });
    psd0.add<bool>("Pointing", true);
    psd0.add<std::vector<double>>("RMax", {120});
    psd0.add<std::vector<double>>("RMin", {60});
    psd0.add<bool>("RandomShoot", false);
    psd0.add<std::vector<double>>("ZMax", {321});
    psd0.add<std::vector<double>>("ZMin", {320});
    psd0.add<bool>("UseDeltaT", false);
    psd0.add<std::vector<double>>("TMin", {0.});
    psd0.add<std::vector<double>>("TMax", {0.05});
    psd0.add<std::vector<double>>("OffsetFirst", {0.});
    desc.add<edm::ParameterSetDescription>("PGunParameters", psd0);
  }
  desc.addUntracked<int>("Verbosity", 0);
  desc.addUntracked<unsigned int>("firstRun", 1);
  desc.add<std::string>("psethack", "random particles in phi and r windows");
  descriptions.add("MultiCloseByParticleGunProducer", desc);
}

void MultiCloseByParticleGunProducer::produce(Event& e, const EventSetup& es) {
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  if (fVerbosity > 0) {
    LogDebug("MultiCloseByParticleGunProducer")
        << " MultiCloseByParticleGunProducer : Begin New Event Generation" << endl;
  }
  fEvt = new HepMC::GenEvent();

  auto const& field = es.getData(m_fieldToken);

  int barcode = 1;
  unsigned int numParticles = fRandomShoot ? CLHEP::RandFlat::shoot(engine, 1, fNParticles) : fNParticles;
  for (size_t ip = 0; ip < fPartIDs.size(); ip++) {
    double phi = CLHEP::RandFlat::shoot(engine, fPhiMin[ip], fPhiMax[ip]);
    double fZ = CLHEP::RandFlat::shoot(engine, fZMin[ip], fZMax[ip]);
    double fR, fEta;
    double fT;

    if (!fControlledByEta) {
      fR = CLHEP::RandFlat::shoot(engine, fRMin[ip], fRMax[ip]);
      fEta = asinh(fZ / fR);
    } else {
      fEta = CLHEP::RandFlat::shoot(engine, fEtaMin[ip], fEtaMax[ip]);
      fR = (fZ / sinh(fEta));
    }

    if (fUseDeltaT) {
      fT = CLHEP::RandFlat::shoot(engine, fTMin[ip], fTMax[ip]);
    } else {
      fT = 0.;
    }

    double tmpPhi = phi;
    double tmpR = fR;

    int PartID = fPartIDs[ip];
    const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID)));
    double mass = PData->mass().value();

    double mom, px, py, pz;
    double energy;
    double fVar = CLHEP::RandFlat::shoot(engine, fVarMin[ip], fVarMax[ip]);
    if (!fFlatPtGeneration) {
      double mom2 = fVar * fVar - mass * mass;
      mom = 0.;
      if (mom2 > 0.) {
        mom = sqrt(mom2);
      }
      px = 0.;
      py = 0.;
      pz = mom;
      energy = fVar;
    } else {
      double theta = 2. * atan(exp(-fEta));
      mom = fVar / sin(theta);
      px = fVar * cos(phi);
      py = fVar * sin(phi);
      pz = mom * cos(theta);
      double energy2 = mom * mom + mass * mass;
      energy = sqrt(energy2);
    }
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

    // compute correct path assuming uniform magnetic field in CMS
    double pathLength = 0.;
    const double speed = p.pz() / p.e() * c_light / CLHEP::cm;
    if (PData->charge()) {
      // Radius [cm] = P[GeV/c] * 10^9 / (c[mm/ns] * 10^6 * q[C] * B[T]) * 100[cm/m]
      const double radius = std::sqrt(p.px() * p.px() + p.py() * p.py()) * std::pow(10, 5) /
                            (c_light * field.inTesla({0.f, 0.f, 0.f}).z());  // cm
      const double arc = 2 * asinf(std::sqrt(x * x + y * y) / (2 * radius)) * radius;
      pathLength = std::sqrt(arc * arc + fZ * fZ);
    } else {
      pathLength = std::sqrt(x * x + y * y + fZ * fZ);
    }

    // if not pointing time doesn't mean a lot, keep the old way
    const double pathTime = fPointing ? (pathLength / speed) : (std::sqrt(x * x + y * y + fZ * fZ) / speed);
    double timeOffset = fOffsetFirst[ip] + (pathTime + ip * fT) * CLHEP::ns * c_light;

    HepMC::GenVertex* Vtx =
        new HepMC::GenVertex(HepMC::FourVector(x * CLHEP::cm, y * CLHEP::cm, fZ * CLHEP::cm, timeOffset));

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
    LogDebug("MultiCloseByParticleGunProducer") << " MultiCloseByParticleGunProducer : Event Generation Done " << endl;
  }
}
