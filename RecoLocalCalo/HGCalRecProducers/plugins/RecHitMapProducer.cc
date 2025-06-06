// user include files
#include <unordered_map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/Common/interface/MultiSpan.h"

class RecHitMapProducer : public edm::global::EDProducer<> {
public:
  RecHitMapProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<HGCRecHitCollection> hits_ee_token_;
  const edm::EDGetTokenT<HGCRecHitCollection> hits_fh_token_;
  const edm::EDGetTokenT<HGCRecHitCollection> hits_bh_token_;
  const edm::EDGetTokenT<MultiCollectionManager<HGCRecHitCollection>> hgcalToken_;
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_eb_token_;
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_hb_token_;
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_ho_token_;
  bool hgcalOnly_;
};

DEFINE_FWK_MODULE(RecHitMapProducer);

using DetIdRecHitMap = std::unordered_map<DetId, const unsigned int>;

RecHitMapProducer::RecHitMapProducer(const edm::ParameterSet& ps)
    : hgcalToken_{consumes<MultiCollectionManager<HGCRecHitCollection>>(
          ps.getParameter<edm::InputTag>("HGCalMultiRecHits"))},
      hits_eb_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("EBInput"))),
      hits_hb_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HBInput"))),
      hits_ho_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HOInput"))),
      hgcalOnly_(ps.getParameter<bool>("hgcalOnly")) {
  produces<DetIdRecHitMap>("hgcalRecHitMap");
  if (!hgcalOnly_)
    produces<DetIdRecHitMap>("barrelRecHitMap");
}

void RecHitMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("HGCalMultiRecHits", {"hgcalRecHitMultiCollectionProducer", ""});
  desc.add<edm::InputTag>("EBInput", {"particleFlowRecHitECAL", ""});
  desc.add<edm::InputTag>("HBInput", {"particleFlowRecHitHBHE", ""});
  desc.add<edm::InputTag>("HOInput", {"particleFlowRecHitHO", ""});
  desc.add<bool>("hgcalOnly", true);
  descriptions.add("recHitMapProducer", desc);
}

void RecHitMapProducer::produce(edm::StreamID, edm::Event& evt, const edm::EventSetup& es) const {
  auto hitMapHGCal = std::make_unique<DetIdRecHitMap>();
  auto const& mgr = evt.get(hgcalToken_);
  auto flat = mgr.makeFlatView();  // by value
  for (unsigned int i = 0; i < flat.size(); ++i) {
    hitMapHGCal->emplace(flat[i].detid(), i);
  }

  // TODO may be worth to avoid dependency on the order
  // of the collections, maybe using a map
  edm::MultiSpan<HGCRecHit> rechitSpan;
  rechitSpan.add(*ee_hits);
  rechitSpan.add(*fh_hits);
  rechitSpan.add(*bh_hits);

  for (unsigned int i = 0; i < rechitSpan.size(); ++i) {
    const auto recHitDetId = rechitSpan[i].detid();
    hitMapHGCal->emplace(recHitDetId, i);
  }

  evt.put(std::move(hitMapHGCal), "hgcalRecHitMap");

  if (!hgcalOnly_) {
    auto hitMapBarrel = std::make_unique<DetIdRecHitMap>();
    edm::MultiSpan<reco::PFRecHit> barrelRechitSpan;
    barrelRechitSpan.add(evt.get(barrel_hits_token_[0]));
    barrelRechitSpan.add(evt.get(barrel_hits_token_[1]));
    for (unsigned int i = 0; i < barrelRechitSpan.size(); ++i) {
      const auto recHitDetId = barrelRechitSpan[i].detId();
      hitMapBarrel->emplace(recHitDetId, i);
    }
    evt.put(std::move(hitMapBarrel), "barrelRecHitMap");
  }
}
