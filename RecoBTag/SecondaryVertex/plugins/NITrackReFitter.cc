// -*- C++ -*-
//
// Package:    TrackReFitter
// Class:      TrackReFitter
// 
/**\class TrackReFitter TrackReFitter.cc Validation/RecoB/src/TrackReFitter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dominik Nowatschin,68/102,2978,
//         Created:  Mon Jan 13 14:35:23 CET 2014
// $Id$
//
//


// system include files
#include <memory>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// #include "TH2.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerAlgorithm.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToBeamLine.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "TrackingTools/TrackFitters/interface/TrajectoryFitterRecord.h" 
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 

#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"


#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TracksClusteringFromDisplacedSeed.h"

// #define TRACKREFIT_DBG 0



// #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <math.h>

namespace ni_analyzer
{
    using namespace reco;
    typedef std::pair<Trajectory*, std::pair<reco::Track*,PropagationDirection> > AlgoProduct; 
    typedef std::vector< AlgoProduct >  AlgoProductCollection;
    typedef edm::RefToBase<TrajectorySeed> SeedRef;

    bool buildTrack (const TrajectoryFitter * theFitter, const Propagator * thePropagator, std::vector<reco::Track*> & algoResults, TransientTrackingRecHit::RecHitContainer& hits, TrajectoryStateOnSurface& theTSOS, const TrajectorySeed& seed, float ndof, const reco::BeamSpot& bs, SeedRef seedRef, int qualityMask,signed char nLoops, reco::TrackBase::TrackAlgorithm algo_, bool geometricInnerState_);

    TrajectoryStateOnSurface getInitialState(const reco::Track * theT, TransientTrackingRecHit::RecHitContainer& hits, const TrackingGeometry * theG, const MagneticField * theMF);

    reco::VertexCompositePtrCandidate refitVertex(reco::VertexCompositePtrCandidate const & vertex, edm::Handle<BeamSpot> beamSpot, std::vector<reco::Track const*> & tracks, edm::ESHandle<TransientTrackBuilder> trackBuilder);
}



class NITrackReFitter : public edm::EDProducer {
public:

    explicit NITrackReFitter(const edm::ParameterSet&);
    ~NITrackReFitter();
    
    virtual void produce(edm::Event&, const edm::EventSetup&);

private:


    edm::EDGetTokenT<edm::View<reco::Candidate> > tokenSecondaryVertexCollection_;
    edm::EDGetTokenT<reco::VertexCollection > token_primaryVertex_;
    edm::EDGetTokenT<reco::BeamSpot> bsSrc_;

    std::vector<TH1D*> SVHistosRho_;
};

NITrackReFitter::NITrackReFitter(const edm::ParameterSet& pSet) 
{
   tokenSecondaryVertexCollection_ = consumes<edm::View<reco::Candidate> >(pSet.getParameter<edm::InputTag>("secondaryVertices"));
   token_primaryVertex_ = consumes<reco::VertexCollection> (pSet.getParameter<edm::InputTag>("primaryVertices"));
   bsSrc_ = consumes<reco::BeamSpot> (pSet.getParameter<edm::InputTag>("beamSpot"));
   produces<std::vector<reco::VertexCompositePtrCandidate> >();
   produces<std::vector<reco::Vertex> >();
//    edm::Service<TFileService> fs;
//    for (unsigned int i = 0; i < secondaryVertexTags_.size(); ++i){
// 	   TH1D *SVhistoRho = fs->make<TH1D>(secondaryVertexTags_[i].label().c_str(), "Secondary Vertices Rho Projection", 300, 0, 30);
// 	   TString title = secondaryVertexTags_[i].label()+";Rho;No. of Events";
// 	   SVhistoRho->SetTitle(title);
// 	   SVHistosRho_.push_back(SVhistoRho);
//    }

}


NITrackReFitter::~NITrackReFitter()
{
}

void NITrackReFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace reco;
    using namespace edm;

    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken(bsSrc_,recoBeamSpotHandle);
    
    edm::Handle<edm::View<reco::Candidate> > secondaryVertices;
    iEvent.getByToken(tokenSecondaryVertexCollection_, secondaryVertices);

    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(token_primaryVertex_,primaryVertices);

    edm::ESHandle<TrackerGeometry> trackingGeometry ;
    iSetup.get<TrackerDigiGeometryRecord>().get(trackingGeometry);

    const reco::Vertex & primVert = *primaryVertices->begin();
    
    
    /// declare the fitters, builders, etc. here; so far, everything is hard-coded, put this in a cfi file later

    edm::ESHandle<MagneticField> theMF;
    edm::ESHandle<TrajectoryFitter> theFitter;
    edm::ESHandle<Propagator> thePropagator;
// 	edm::ESHandle<MeasurementTracker>  theMeasTk;
    edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
// 	getFromES(setup,theG,theMF,theFitter,thePropagator,theMeasTk,theBuilder);
    
// 	if (conf_.exists("useSimpleMF")) useSimpleMF = conf_.getParameter<bool>("useSimpleMF");
    std::string mfName = "";
// 	if (useSimpleMF){
// 		mfName = conf_.getParameter<std::string>("SimpleMagneticField"); 
// 	}
    iSetup.get<IdealMagneticFieldRecord>().get(mfName, theMF);

    std::string fitterName = "KFFittingSmootherWithOutliersRejectionAndRK";			// = conf_.getParameter<std::string>("Fitter");   
    iSetup.get<TrajectoryFitter::Record>().get(fitterName,theFitter);

    std::string propagatorName = "RungeKuttaTrackerPropagator";			// = conf_.getParameter<std::string>("Propagator");   
    iSetup.get<TrackingComponentsRecord>().get(propagatorName,thePropagator);

// 	edm::ESHandle<NavigationSchool> theSchool;
// 	std::string theNavigationSchool ="SimpleNavigationSchool";
// 	if (conf_.exists("NavigationSchool")) theNavigationSchool= conf_.getParameter<std::string>("NavigationSchool");
// 	else edm::LogWarning("TrackProducerBase")<<" NavigationSchool parameter not set. secondary hit pattern will not be filled.";
// 	if (theNavigationSchool!=""){
// 	iSetup.get<NavigationSchoolRecord>().get(theNavigationSchool, theSchool);

// 	std::string measTkName = "";			// = conf_.getParameter<std::string>("MeasurementTracker");
// 	iSetup.get<CkfComponentsRecord>().get(measTkName,theMeasTk);

    std::string builderName = "WithAngleAndTemplate";			// = conf_.getParameter<std::string>("TTRHBuilder");   
    iSetup.get<TransientRecHitRecord>().get(builderName,theBuilder);
    
    //    std::auto_ptr<std::vector<reco::VertexCompositePtrCandidate> > newVertices (new std::vector<reco::VertexCompositePtrCandidate>);
    auto newVertices = std::make_unique<std::vector<reco::VertexCompositePtrCandidate> >();
    auto newSVs = std::make_unique<std::vector<reco::Vertex> >();

    /// loop over all vertices selected for re-fitting (normally the ones close to pixel layers)
    
    int count = 0;
    for (edm::View<reco::Candidate>::const_iterator iVertex = secondaryVertices->begin(); iVertex != secondaryVertices->end(); ++iVertex, ++count){
    
        GlobalPoint sv(iVertex->vx(), iVertex->vy(), iVertex->vz());
        GlobalPoint pv(primVert.position().x(), primVert.position().y(), primVert.position().z());

#ifdef TRACKREFIT_DBG
        std::cout << "  Vertex " << count << " at rho / phi / z: " << sv.perp() << " / " << sv.phi() << " / " << sv.z() << std::endl;
#endif


        std::vector<reco::Track*> algoResults;
        std::vector<reco::Track const*> newTracks;  /// FIXME: try to change this to something like vector<PFCandidate> so later you won't get runtime errors/warnings and are able to return a pointer in ni_analyzer::refitVertex
        bool reFit = false;

        /// loop over all daughter tracks of a vertex and check whether they should be re-fitted
        
        for (size_t iTrack = 0; iTrack < iVertex->numberOfSourceCandidatePtrs(); ++iTrack){
#ifdef TRACKREFIT_DBG
            std::cout << "    Track " << iTrack << ":" << std::endl;
#endif
            CandidatePtr trackCand = iVertex->sourceCandidatePtr(iTrack);
            reco::Track const & recoTrack = *trackCand->bestTrack();

#ifdef TRACKREFIT_DBG
            std::cout << "     Innermost hit rho / phi / x / y / z: " << recoTrack.innerPosition().rho() << " / " << recoTrack.innerPosition().phi() << " / " << recoTrack.innerPosition().x() << " / " << recoTrack.innerPosition().y() << " / " << recoTrack.innerPosition().z() << std::endl << std::endl;
#endif

            /// if the innermost hit of a track is closer to the beam axis than the vertex this track is coming from, don't re-fit assuming that this track doesn't come from this vertex anyway
            if (sv.perp() - recoTrack.innerPosition().rho() < 0.2) {

                const TkTransientTrackingRecHitBuilder * builder = dynamic_cast<TkTransientTrackingRecHitBuilder const *>(theBuilder.product());
                assert(builder);
                try{
                    float ndof=0;
                    PropagationDirection seedDir = recoTrack.seedDirection();

                    // WARNING: here we assume that the hits are correcly sorted according to seedDir
                    TransientTrackingRecHit::RecHitContainer hits;
                    
                    /// in the following, the crucial part is '... iHit = recoTrack.recHitsBegin()+1' so we use all hits from this track EXCEPT the innermost one
                    for (trackingRecHit_iterator iHit=recoTrack.recHitsBegin()+1; iHit!=recoTrack.recHitsEnd(); iHit++){
                        if ((**iHit).geographicalId()!=0U)  hits.push_back( (**iHit).cloneForFit(*builder->geometry()->idToDet( (**iHit).geographicalId() ) ) );
                    }
	
                    TrajectoryStateOnSurface theInitialStateForRefitting = ni_analyzer::getInitialState(&recoTrack,hits,trackingGeometry.product(),theMF.product());

                    // the seed has dummy state and hits.What matters for the fitting is the seedDirection;
                    if (seedDir==anyDirection){//if anyDirection the seed direction is not stored in the root file: keep same order
                        throw cms::Exception("TrackProducer") << "ERROR: trying to refit a track which doesn't have a properly filled 'seed direction' data member" << std::endl;
                    }

                    const TrajectorySeed seed = TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), seedDir);
                    reco::TrackBase::TrackAlgorithm algo_=recoTrack.algo();

                    //=====  the hits are in the same order as they were in the track::extra.
                    FitterCloner fc(theFitter.product(),builder);      

                    bool geometricInnerState_ = false;

                    bool ok = ni_analyzer::buildTrack(fc.fitter.get(),thePropagator.product(),algoResults, hits, theInitialStateForRefitting, seed, ndof, *recoBeamSpotHandle, recoTrack.seedRef(),recoTrack.qualityMask(),recoTrack.nLoops(), algo_, geometricInnerState_);
                    if(ok) {
                        reco::Track const * refittedTrack = algoResults.back();
                        ///FIXME: build new PFCandidate from refitted Track here and push back
                        newTracks.push_back(refittedTrack);
                        reFit = true;
//                         PXBDetId originalPxbDetId((*recoTrack.recHitsBegin())->geographicalId());
//                         PXBDetId refitPxbDetId((*refittedTrack->recHitsBegin())->geographicalId());
//                         std::cout << "TEST0" << std::endl;
#ifdef TRACKREFIT_DBG
                        std::cout << "     Track Refit successful!" << std::endl
                            << "       Original track chi2/normalizedChi2: " << recoTrack.chi2() << "/" << recoTrack.normalizedChi2()
                            << " | refitted track chi2/normalizedChi2: " << " | " << refittedTrack->chi2()  << "/" << refittedTrack->normalizedChi2() << std::endl
                            << "       Original track pt/ptError: " << recoTrack.pt() << "/" << recoTrack.ptError()
                            << " | refitted track pt/pterror: " << " | " << refittedTrack->pt() << "/" << refittedTrack->ptError() << std::endl
                            << "       Original track phi/phiError: " << recoTrack.phi() << "/" << recoTrack.phiError()
                            << " | refitted track phi/phierror: " << " | " << refittedTrack->phi() << "/" << refittedTrack->phiError() << std::endl
                            << "       Original track eta/etaError: " << recoTrack.eta() << "/" << recoTrack.etaError()
                            << " | refitted track eta/etaerror: " << " | " << refittedTrack->eta() << "/" << refittedTrack->etaError() << std::endl
                            << "       Original track dxy/dxyError: " << recoTrack.dxy(*recoBeamSpotHandle) << "/" << recoTrack.dxyError()
                            << " | refitted track dxy/dxyerror: " << " | " << refittedTrack->dxy(*recoBeamSpotHandle) << "/" << refittedTrack->dxyError() << std::endl
                            << "       Original track dz/dzError: " << recoTrack.dz(recoBeamSpotHandle->position(recoTrack.vz())) << "/" << recoTrack.dzError()
                            << " | refitted track dz/dzerror: " << " | " << refittedTrack->dz(recoBeamSpotHandle->position(refittedTrack->vz())) << "/" << refittedTrack->dzError() << std::endl;
#endif
                        continue;
//                         std::cout << "     Original track innermost hit position/id/layer | refitted track innermost hit position/id/layer: " << recoTrack.innerPosition().rho() << "/" << originalPxbDetId.subdetId() << "/" << originalPxbDetId.layer() << " | " << refittedTrack->innerPosition().rho() << "/" << refitPxbDetId.subdetId() << "/" << refitPxbDetId.layer() << std::endl;
                    }
#ifdef TRACKREFIT_DBG
                    else std::cout << "     NOT successful!" << std::endl;
#endif
                }catch ( cms::Exception & e){
                    std::cout << "Genexception1: " << e.explainSelf() <<"\n";      
                    throw;
                }
                
            }

#ifdef TRACKREFIT_DBG
            std::cout << "     No tracks refitted, add original Track!" << std::endl;
#endif

            newTracks.push_back(trackCand->bestTrack());
  
        }

        reco::VertexCompositePtrCandidate const * originalVertex = dynamic_cast<reco::VertexCompositePtrCandidate const *>(&*iVertex);
        
        if (reFit) {
#ifdef TRACKREFIT_DBG
            std::cout << "   Refit Vertex!" << std::endl;
#endif

            reco::VertexCompositePtrCandidate newVertex = ni_analyzer::refitVertex(*originalVertex, recoBeamSpotHandle, newTracks, builder);
            
            newVertices->push_back(newVertex);
	    reco::Vertex newSV(newVertex.vertex(), newVertex.error4D(), newVertex.t(), newVertex.vertexChi2(), newVertex.vertexNdof(), 10);
	    newSVs->push_back(newSV);

#ifdef TRACKREFIT_DBG
            if (!newVertex.vertexChi2()) std::cout << "     Vertex fitting NOT successful!" << std::endl;
            else std::cout << "   Vertex fitting successful!" << std::endl
                << "     Original vertex rho: " << originalVertex->vertex().rho()
                << " | refitted vertex rho: " << " | " << newVertex.vertex().rho() << std::endl
                << "     Original vertex chi2/normalizedChi2: " << originalVertex->vertexChi2() << "/" << originalVertex->vertexNormalizedChi2()
                << " | refitted vertex chi2/normalizedChi2: " << " | " << newVertex.vertexChi2()  << "/" << newVertex.vertexNormalizedChi2() << std::endl
                << "     Original vertex x/xCovEl: " << originalVertex->vertex().x() << "/" << originalVertex->vertexCovariance(0,0)
                << " | refitted vertex x/xCovEl: " << " | " << newVertex.vertex().x() << "/" << newVertex.vertexCovariance(0,0) << std::endl
                << "     Original vertex y/yCovEl: " << originalVertex->vertex().y() << "/" << originalVertex->vertexCovariance(1,1)
                << " | refitted vertex y/yCovEl: " << " | " << newVertex.vertex().y() << "/" << newVertex.vertexCovariance(1,1) << std::endl
                << "     Original vertex z/zCovEl: " << originalVertex->vertex().z() << "/" << originalVertex->vertexCovariance(2,2)
                << " | refitted vertex z/zCovEl: " << " | " << newVertex.vertex().z() << "/" << newVertex.vertexCovariance(2,2) << std::endl << std::endl;
// 				<< "     Original vertex z/dxyError: " << originalVertex->dxy(*recoBeamSpotHandle) << "/" << originalVertex->dxyError()
// 				<< " | refitted vertex dxy/dxyerror: " << " | " << newVertex.dxy(*recoBeamSpotHandle) << "/" << newVertex.dxyError() << std::endl
// 				<< "     Original vertex dz/dzError: " << originalVertex->dz(recoBeamSpotHandle->position(originalVertex->vz())) << "/" << originalVertex->dzError()
// 				<< " | refitted vertex dz/dzerror: " << " | " << newVertex.dz(recoBeamSpotHandle.position(newVertex.vz())) << "/" << newVertex.dzError() << std::endl;
#endif
                
        }
        else newVertices->push_back(*originalVertex);

    }
    iEvent.put(std::move(newVertices));
   
   
}

bool ni_analyzer::buildTrack (const TrajectoryFitter * theFitter, const Propagator * thePropagator, std::vector<reco::Track*> & algoResults, TransientTrackingRecHit::RecHitContainer& hits, TrajectoryStateOnSurface& theTSOS, const TrajectorySeed& seed, float ndof, const reco::BeamSpot& bs, SeedRef seedRef, int qualityMask,signed char nLoops, reco::TrackBase::TrackAlgorithm algo_, bool geometricInnerState_) 
{
//   PropagationDirection seedDir = seed.direction();
      
    //perform the fit: the result's size is 1 if it succeded, 0 if fails
    Trajectory && trajTmp = theFitter->fitOne(seed, hits, theTSOS,(nLoops>0) ? TrajectoryFitter::looper : TrajectoryFitter::standard);
    if UNLIKELY(!trajTmp.isValid()) {
        return false;
    }
  
  
    auto theTraj = new Trajectory(std::move(trajTmp));
    theTraj->setSeedRef(seedRef);
  
    ndof = 0;
    for (auto const & tm : theTraj->measurements()) {
        auto const & h = tm.recHitR();
        if (h.isValid()) ndof = ndof + float(h.dimension())*h.weight();  // two virtual calls!
    }
  
    ndof -= 5.f;
    if UNLIKELY(std::abs(theTSOS.magneticField()->nominalValue())<DBL_MIN) ++ndof;  // same as -4
 
    //if geometricInnerState_ is false the state for projection to beam line is the state attached to the first hit: to be used for loopers
    //if geometricInnerState_ is true the state for projection to beam line is the one from the (geometrically) closest measurement to the beam line: to be sued for non-collision tracks
    //the two shouuld give the same result for collision tracks that are NOT loopers
    TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface;
    if (geometricInnerState_) {
        stateForProjectionToBeamLineOnSurface = theTraj->closestMeasurement(GlobalPoint(bs.x0(),bs.y0(),bs.z0())).updatedState();
    } else {
        if (theTraj->direction() == alongMomentum) {
            stateForProjectionToBeamLineOnSurface = theTraj->firstMeasurement().updatedState();
        } else { 
            stateForProjectionToBeamLineOnSurface = theTraj->lastMeasurement().updatedState();
        }
    }

    if UNLIKELY(!stateForProjectionToBeamLineOnSurface.isValid()){
//     edm::LogError("CannotPropagateToBeamLine")<<"the state on the closest measurement isnot valid. skipping track.";
        delete theTraj;
        return false;
    }
    const FreeTrajectoryState & stateForProjectionToBeamLine=*stateForProjectionToBeamLineOnSurface.freeState();
  
//   LogDebug("TrackProducer") << "stateForProjectionToBeamLine=" << stateForProjectionToBeamLine;
  
    TSCBLBuilderNoMaterial tscblBuilder;
    TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine,bs);
  
    if UNLIKELY(!tscbl.isValid()) {
        delete theTraj;
        return false;
    }
  
    GlobalPoint v = tscbl.trackStateAtPCA().position();
    math::XYZPoint  pos( v.x(), v.y(), v.z() );
    GlobalVector p = tscbl.trackStateAtPCA().momentum();
    math::XYZVector mom( p.x(), p.y(), p.z() );
  
//   LogDebug("TrackProducer") << "pos=" << v << " mom=" << p << " pt=" << p.perp() << " mag=" << p.mag();
  
    std::auto_ptr<reco::Track> theTrack(new reco::Track(theTraj->chiSquared(), int(ndof), pos, mom, tscbl.trackStateAtPCA().charge(), tscbl.trackStateAtPCA().curvilinearError(), algo_));
  
    theTrack->setQualityMask(qualityMask);
    theTrack->setNLoops(nLoops);
    
    algoResults.push_back(theTrack.get());
  

    return true;
} 


TrajectoryStateOnSurface ni_analyzer::getInitialState(const reco::Track * theT, TransientTrackingRecHit::RecHitContainer& hits, const TrackingGeometry * theG, const MagneticField * theMF){

    TrajectoryStateOnSurface theInitialStateForRefitting;
    //the starting state is the state closest to the first hit along seedDirection.
    //avoiding to use transientTrack, it should be faster;
    TrajectoryStateOnSurface innerStateFromTrack=trajectoryStateTransform::innerStateOnSurface(*theT,*theG,theMF);
    TrajectoryStateOnSurface outerStateFromTrack=trajectoryStateTransform::outerStateOnSurface(*theT,*theG,theMF);
    TrajectoryStateOnSurface initialStateFromTrack = 
        ( (innerStateFromTrack.globalPosition()-hits.front()->det()->position()).mag2() <
        (outerStateFromTrack.globalPosition()-hits.front()->det()->position()).mag2() ) ? 
        innerStateFromTrack: outerStateFromTrack;       
  
    // error is rescaled, but correlation are kept.
    initialStateFromTrack.rescaleError(100);
    return initialStateFromTrack;
}



reco::VertexCompositePtrCandidate ni_analyzer::refitVertex(reco::VertexCompositePtrCandidate const & vertex, edm::Handle<BeamSpot> beamSpot, std::vector<reco::Track const*> & tracks, edm::ESHandle<TransientTrackBuilder> trackBuilder)
{
#ifdef TRACKREFIT_DBG
            std::cout << "   Refit Vertex start" << std::endl;
#endif

    GlobalPoint originalVertex(vertex.vx(), vertex.vy(), vertex.vz());
    double sigmacut = 3.0;
    double Tini     = 256.;
    double ratio    = 0.25;
// 	VertexDistance3D vdist;
// 	VertexDistanceXY vdist2d;
// 	MultiVertexFitter theMultiVertexFitter;
    AdaptiveVertexFitter theAdaptiveFitter(GeometricAnnealing(sigmacut, Tini, ratio), DefaultLinearizationPointFinder(), KalmanVertexUpdator<5>(), KalmanVertexTrackCompatibilityEstimator<5>(), KalmanVertexSmoother());

#ifdef TRACKREFIT_DBG
            std::cout << "   theAdaptiveFitter" << std::endl;
#endif

        
    std::vector<TransientTrack> tts;
    
    //Fill transient track vector 
    for(size_t iTrack = 0; iTrack != tracks.size(); ++iTrack) {
        TransientTrack tt = trackBuilder->build(tracks[iTrack]);
        tt.setBeamSpot(*beamSpot);
        tts.push_back(tt);
    }
#ifdef TRACKREFIT_DBG
            std::cout << "   tts filled" << std::endl;
#endif


    TransientVertex singleFitVertex = theAdaptiveFitter.vertex(tts,originalVertex);

    /// originally wanted to have 'std::auto_ptr<reco::VertexCompositePtrCandidate> refittedVertex', so one can easily check whether vertex was fitted by checking if this function returns a null pointer; however, this throws a runtime error (see log from 24/11/14) which is actually still present in the current setup
    
    reco::VertexCompositePtrCandidate refittedVertex;

#ifdef TRACKREFIT_DBG
    std::cout << "   vertex single fit done" << std::endl;
#endif


    if(singleFitVertex.isValid()){
#ifdef TRACKREFIT_DBG
    std::cout << "   single vertex fit valid" << std::endl;
#endif

		VertexCompositePtrCandidate dummyVertex(singleFitVertex);
#ifdef TRACKREFIT_DBG
    std::cout << "   dummy vertex" << std::endl;
#endif

// 		std::auto_ptr<reco::VertexCompositePtrCandidate> refittedVertex(new reco::VertexCompositePtrCandidate(singleFitVertex));
		refittedVertex = dummyVertex;
#ifdef TRACKREFIT_DBG
    std::cout << "   refitted vertex" << std::endl;
#endif

// 		*refittedVertex = dummyVertex;
// 		return refittedVertex;
	}

#ifdef TRACKREFIT_DBG
            std::cout << "   Refitted Vertex" << std::endl;
#endif

	
	return refittedVertex;
        
}


DEFINE_FWK_MODULE(NITrackReFitter);
