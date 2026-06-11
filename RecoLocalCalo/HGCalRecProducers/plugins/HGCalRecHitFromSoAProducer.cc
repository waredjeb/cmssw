#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

class HGCalRecHitFromSoAProducer : public edm::stream::EDProducer<> {
public:
  HGCalRecHitFromSoAProducer(edm::ParameterSet const& config)
      : getTokenSoARecHits_(consumes(config.getParameter<edm::InputTag>("src"))) {
    produces<HGCRecHitCollection>();
  }

  ~HGCalRecHitFromSoAProducer() override = default;

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override {
    auto const& soaRecHits = iEvent.get(getTokenSoARecHits_);
    auto const soaView = soaRecHits.view();

    auto output = std::make_unique<HGCRecHitCollection>();
    output->reserve(soaRecHits->metadata().size());

    for (int i = 0; i < soaRecHits->metadata().size(); ++i) {
      DetId detid(soaView[i].detid());

      HGCRecHit rechit(detid, soaView[i].energy(), soaView[i].time(), 0, 0);
      rechit.setTimeError(soaView[i].timeError());

      output->push_back(rechit);
    }

    iEvent.put(std::move(output));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hgcalRecHits"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  edm::EDGetTokenT<HGCalSoARecHitsHostCollection> const getTokenSoARecHits_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalRecHitFromSoAProducer);
