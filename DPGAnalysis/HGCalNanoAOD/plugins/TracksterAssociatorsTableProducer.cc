#include "PhysicsTools/NanoAOD/interface/AssociationMapFlatTableProducer.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

typedef AssociationOneToOneFlatTableProducer<TICLAssociationMapOneToOneFraction<SimCluster, CaloParticle>>
    SimClusterCaloParticleFractionFlatTableProducer;

typedef AssociationOneToManyFlatTableProducer<TICLAssociationMapOneToManySharedEnergyScore<ticl::Trackster, ticl::Trackster>>
    TracksterTracksterEnergyScoreFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimClusterCaloParticleFractionFlatTableProducer);
DEFINE_FWK_MODULE(TracksterTracksterEnergyScoreFlatTableProducer);
