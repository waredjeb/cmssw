#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
typedef SimpleCollectionFlatTableProducer<HGCRecHit> RecHitsCollectionTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecHitsCollectionTableProducer);
