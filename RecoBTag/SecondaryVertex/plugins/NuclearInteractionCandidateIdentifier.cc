#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include <set>





class NuclearInteractionCandidateIdentifier : public edm::EDProducer {
    public:
	NuclearInteractionCandidateIdentifier(const edm::ParameterSet &params);


	virtual void produce(edm::Event &event, const edm::EventSetup &es);

    private:

	edm::EDGetTokenT<reco::VertexCollection> tokenPrimaryVertexCollection;
	edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate> > tokenSecondaryVertexCollection;
// 	edm::InputTag beamSpotCollection;
	edm::ParameterSet selectionCriteria;
};

NuclearInteractionCandidateIdentifier::NuclearInteractionCandidateIdentifier(const edm::ParameterSet &params) :
// 	beamSpotCollection           (params.getParameter<edm::InputTag>("beamSpot"))
	selectionCriteria(params.getParameter<edm::ParameterSet>("selection"))
{
	tokenPrimaryVertexCollection = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
	tokenSecondaryVertexCollection = consumes<edm::View<reco::VertexCompositePtrCandidate> >(params.getParameter<edm::InputTag>("secondaryVertices"));
	produces<edm::PtrVector<reco::Candidate> >();
	//produces<std::vector<reco::VertexCompositePtrCandidate> >();
}


void NuclearInteractionCandidateIdentifier::produce(edm::Event &event, const edm::EventSetup &es)
{
// 	using namespace reco;

	typedef reco::Candidate::Point Point;
	typedef reco::Candidate::Vector Vector;

	edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > secondaryVertices;
	event.getByToken(tokenSecondaryVertexCollection, secondaryVertices);
	
	edm::Handle<reco::VertexCollection> primaryVertices;
	event.getByToken(tokenPrimaryVertexCollection, primaryVertices);

	//	std::auto_ptr<edm::PtrVector<reco::Candidate> > niVertices(new edm::PtrVector<reco::Candidate>);
	auto niVertices = std::make_unique<edm::PtrVector<reco::Candidate> >();
	//auto recoVertices = std::make_unique<std::vector<reco::VertexCompositePtrCandidate> >();
	//	std::auto_ptr<std::vector<reco::VertexCompositePtrCandidate> > recoVertices(new std::vector<reco::VertexCompositePtrCandidate>);
	// 		const reco::Vertex &pv = (*primaryVertices)[0];
	
	// 		edm::Handle<BeamSpot> beamSpot;
	// 		event.getByLabel(beamSpotCollection, beamSpot);


	for(unsigned int ivtx=0; ivtx < secondaryVertices->size(); ivtx++){
// 		std::cout << std::endl << "======NEW VERTEX========" << std::endl << std::endl;
		const reco::VertexCompositePtrCandidate & sv = (*secondaryVertices)[ivtx];
		// 			GlobalPoint pvPos(pv.position().x(),pv.position().y(),pv.position().z());
		GlobalPoint svPos(sv.vx(),sv.vy(),sv.vz());
		
// 		if (svPos.rho() != gsv.perp()) std::cout << "Rho values are not equal, svPos.rho()-gsv.perp() = " << svPos.rho()-gsv.perp() << std::endl;
		float mass=sv.mass();
		
		bool isNI = false;
		
		// check whether vertex is on detector material
		
		if (selectionCriteria.exists("position")){
			std::vector<double> layerCuts = selectionCriteria.getParameter<std::vector<double> >("position");
			
			for (unsigned int iLayer = 0; iLayer < layerCuts.size(); ++iLayer){
				if (iLayer % 2 == 0 && svPos.perp() >= layerCuts[iLayer]) isNI = true;
				if (iLayer % 2 != 0 && svPos.perp() > layerCuts[iLayer]) isNI = false;					
			}
			
// 			if (isNI)
// 			std::cout << "NI identified at rho = " << svPos.rho() << " in event " << event.id() << std::endl;
		}
		
		// if it is close to detector material, check for other criteria

		if (selectionCriteria.exists("maxZ") && isNI){
			float z = std::abs(svPos.z());
			double maxZ = selectionCriteria.getParameter<double>("maxZ");
			if (z > maxZ) isNI = false;
		}
		
		if (selectionCriteria.exists("minNctau") && isNI){
			const reco::Vertex &pv = (*primaryVertices)[0];
			GlobalPoint pvPos(pv.position().x(),pv.position().y(),pv.position().z());
			float pt=sv.pt();
			float gamma=pt/mass;
			GlobalVector flightDir = svPos-pvPos;
			float flightDistance2D = flightDir.perp();
			float Bctau = 0.05;		// c*tau for B hadron
			
			float nctau = flightDistance2D/(gamma*Bctau);	// number of c*taus for potential B hadron
			double minNctau = selectionCriteria.getParameter<double>("minNctau");
			if (nctau < minNctau) isNI = false;
		}
		
		if (selectionCriteria.exists("minMass") && isNI){
			double minMass = selectionCriteria.getParameter<double>("minMass");
			if (mass < minMass) isNI = false;
		}
		
		if (selectionCriteria.exists("maxMass") && isNI){
			double maxMass = selectionCriteria.getParameter<double>("maxMass");
			if (mass > maxMass) isNI = false;
		}
		
		if (selectionCriteria.exists("minNtracks") && isNI){
			int ntracks = sv.numberOfDaughters();
			int minNtracks = selectionCriteria.getParameter<int>("minNtracks");
			if (ntracks < minNtracks) isNI = false;
		}
		
		if (selectionCriteria.exists("maxNtracks") && isNI){
			int ntracks = sv.numberOfDaughters();
			int maxNtracks = selectionCriteria.getParameter<int>("maxNtracks");
			if (ntracks > maxNtracks) isNI = false;	
		}
		
		if (selectionCriteria.exists("minTrack3DipSig") && isNI){
			
			double minTrack3DipSig = selectionCriteria.getParameter<double>("minTrack3DipSig");
			
			edm::ESHandle<TransientTrackBuilder> builder;
			es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
			
			const reco::Vertex &primVert = (*primaryVertices)[0];
			
			GlobalPoint svPos(sv.vx(),sv.vy(),sv.vz());
			GlobalPoint pvPos(primVert.position().x(), primVert.position().y(), primVert.position().z());
			GlobalVector flightDirection = svPos-pvPos;
			
			std::set<double> track3Dips;
			for (size_t iDaughter = 0; iDaughter != sv.numberOfSourceCandidatePtrs(); ++iDaughter){
				
				reco::CandidatePtr trackCand = sv.sourceCandidatePtr(iDaughter);
				reco::TransientTrack transientTrack = builder->build(trackCand);
			
				Measurement1D ip3d = IPTools::signedImpactParameter3D(transientTrack, flightDirection, primVert).second;
// 				Measurement1D ip2d = IPTools::signedTransverseImpactParameter(transientTrack, flightDirection, primVert).second;
				track3Dips.insert(ip3d.significance());
			}
			
// 			std::cout << "Vertex track3Dips: ";
// 			for (std::set<double>::const_iterator i = track3Dips.begin(); i != track3Dips.end(); ++i)
// 				std::cout << *i << " ";
			
// 			std::cout << std::endl;
			
			double highestTrackSig = *track3Dips.end();
			
			if (highestTrackSig < minTrack3DipSig) isNI = false;
		}
		
		if (selectionCriteria.exists("maxTrack3DipSig") && isNI){
			
			double maxTrack3DipSig = selectionCriteria.getParameter<double>("maxTrack3DipSig");
			
			edm::ESHandle<TransientTrackBuilder> builder;
			es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
			
			const reco::Vertex &primVert = (*primaryVertices)[0];
			
			GlobalPoint svPos(sv.vx(),sv.vy(),sv.vz());
			GlobalPoint pvPos(primVert.position().x(), primVert.position().y(), primVert.position().z());
			GlobalVector flightDirection = svPos-pvPos;
			
			std::set<double> track3Dips;
			for (size_t iDaughter = 0; iDaughter != sv.numberOfSourceCandidatePtrs(); ++iDaughter){
				
				reco::CandidatePtr trackCand = sv.sourceCandidatePtr(iDaughter);
				reco::TransientTrack transientTrack = builder->build(trackCand);
			
				Measurement1D ip3d = IPTools::signedImpactParameter3D(transientTrack, flightDirection, primVert).second;
// 				Measurement1D ip2d = IPTools::signedTransverseImpactParameter(transientTrack, flightDirection, primVert).second;
				track3Dips.insert(ip3d.significance());
			}
			
// 			std::cout << "Vertex track3Dips: ";
// 			for (std::set<double>::const_iterator i = track3Dips.begin(); i != track3Dips.end(); ++i)
// 				std::cout << *i << " ";
			
// 			std::cout << std::endl;
			
			double highestTrackSig = *track3Dips.end();
			
			if (highestTrackSig > maxTrack3DipSig) isNI = false;
		}	

		
		if(isNI) {
		  niVertices->push_back(secondaryVertices->ptrAt(ivtx));
		  //		  recoVertices->push_back(sv);
		}
		
		
		
	}
	
	if (selectionCriteria.exists("distToNI")){
		double distToNI = selectionCriteria.getParameter<double>("distToNI");
		for (size_t ivtx = 0; ivtx < secondaryVertices->size(); ++ivtx){
			for (size_t iNI = 0; iNI < niVertices->size(); ++iNI) {
				VertexDistance3D dist;
				reco::VertexCompositePtrCandidate const & sv = (*secondaryVertices)[ivtx];
				reco::CandidatePtr ni = (*niVertices)[iNI];
				VertexState svs(RecoVertex::convertPos(sv.vertex()),RecoVertex::convertError(sv.vertexCovariance()));
				VertexState nis(RecoVertex::convertPos(ni->vertex()),RecoVertex::convertError(ni->vertexCovariance()));
				Measurement1D d = dist.distance(svs, nis);
				if( d.significance() < 3/distToNI and d.value() < distToNI )
				{
				  niVertices->push_back(secondaryVertices->ptrAt(ivtx));
					break;
				}
			}
		}
	}
	
	
	//	event.put(niVertices);
	event.put(std::move(niVertices));
	//	event.put(std::move(recoVertices));
}

DEFINE_FWK_MODULE(NuclearInteractionCandidateIdentifier);
