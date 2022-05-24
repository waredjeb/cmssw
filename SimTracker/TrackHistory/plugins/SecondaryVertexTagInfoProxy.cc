
// system include files
#include <memory>

// user include files
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class SecondaryVertexTagInfoProxy : public edm::global::EDProducer<>
{
public:

    explicit SecondaryVertexTagInfoProxy(const edm::ParameterSet&);
    typedef edm::AssociationMap<edm::OneToMany<reco::SecondaryVertexTagInfoCollection, reco::VertexCollection> > AssocSVTagInfoVertex;

private:

    void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

    edm::EDGetTokenT<reco::SecondaryVertexTagInfoCollection> svTagInfoCollection_;
};


SecondaryVertexTagInfoProxy::SecondaryVertexTagInfoProxy(const edm::ParameterSet& config)
{
    // Get the cfg parameter
    svTagInfoCollection_ = consumes<reco::SecondaryVertexTagInfoCollection>(config.getUntrackedParameter<edm::InputTag> ( "svTagInfoProducer" ));

    // Declare the type of objects to be produced.
    produces<reco::VertexCollection>();
    produces<AssocSVTagInfoVertex >();
}


void SecondaryVertexTagInfoProxy::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const
{
    // Vertex collection
    edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfoCollection;
    event.getByToken(svTagInfoCollection_, svTagInfoCollection);

    // Auto pointers to the collection to be added to the event
    std::unique_ptr<reco::VertexCollection> proxy (new reco::VertexCollection);
    auto assoc = std::make_unique<AssocSVTagInfoVertex>(&event.productGetter());

    // Get a reference before to put in the event
    reco::VertexRefProd vertexRefProd = event.getRefBeforePut<reco::VertexCollection>();

    // General index
    std::size_t index = 0;

    // Loop over SecondaryVertexTagInfo collection
    for (std::size_t svIndex = 0; svIndex < svTagInfoCollection->size(); ++svIndex)
    {
        // Reference to svTagInfo
        reco::SecondaryVertexTagInfoRef svTagInfo(svTagInfoCollection, svIndex);

        // Loop over the vertexes and add them to the new collection
        for (unsigned int vIndex = 0; vIndex < svTagInfo->nVertices(); ++vIndex)
        {
            proxy->push_back(svTagInfo->secondaryVertex(vIndex));
	    assoc->insert(edm::Ref<std::vector<reco::SecondaryVertexTagInfo> >(svTagInfoCollection, svIndex),
			  edm::Ref<std::vector<reco::Vertex> >(vertexRefProd, index));
            ++index;
        }
    }

    // Adding the collection to the event
    event.put(std::move(proxy));
    event.put(std::move(assoc));
}


//define this as a plug-in
DEFINE_FWK_MODULE(SecondaryVertexTagInfoProxy);
