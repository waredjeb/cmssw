#include "FWCore/Framework/interface/ESHandle.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/Common/interface/MultiSpan.h"

class RecHitsExtraTableProducer : public edm::stream::EDProducer<> {
public:
  RecHitsExtraTableProducer(edm::ParameterSet const& params)
      : skipNonExistingSrc_(params.getParameter<bool>("skipNonExistingSrc")),
        tableName_(params.getParameter<std::string>("tableName")),
        precision_(params.getParameter<int>("precision")),
        geom_token_(esConsumes())
  {
  for (auto const &tag : params.getParameter<std::vector<edm::InputTag>>("recHits")) {
    hits_token_.emplace_back(consumes<std::vector<HGCRecHit>>(tag));
  }
    produces<nanoaod::FlatTable>(tableName_);
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& es) override {
    //Layer Clusters time value map
    const auto& geom = es.getData(geom_token_);
    rhtools_.setGeometry(geom);

    static constexpr float default_value = std::numeric_limits<float>::quiet_NaN();
    std::vector<edm::Handle<std::vector<HGCRecHit>>> rechits_h(hits_token_.size());
    edm::MultiSpan<HGCRecHit> recHitsManager;
    for (unsigned int i = 0; i < hits_token_.size(); ++i) {
      iEvent.getByToken(hits_token_[i], rechits_h[i]);
      //Fill MultiSpan
      recHitsManager.add(*rechits_h[i]);
    }
    auto const nRecHits = recHitsManager.size();

    std::vector<float> posX(nRecHits, default_value);
    std::vector<float> posY(nRecHits, default_value);
    std::vector<float> posZ(nRecHits, default_value);
    std::vector<float> energy(nRecHits, default_value);
    std::vector<int>   thickness(nRecHits, 0);
    std::vector<float> layer(nRecHits, default_value);
    std::vector<bool>  isHalfCell(nRecHits,  false);

    // initialize to quiet Nans
    if (!(skipNonExistingSrc_)) {
      for (size_t i = 0; i < recHitsManager.size(); ++i) {
        auto const& rechit = recHitsManager[i];
        auto const& detid = rechit.detid();
        const GlobalPoint position(rhtools_.getPosition(rechit.detid()));
        posX[i] = position.x();
        posY[i] = position.y();
        posZ[i] = position.z();
        energy[i] = rechit.energy();
        thickness[i] = rhtools_.getSiThickness(detid); //0 if scintillator 
        layer[i] = rhtools_.getLayerWithOffset(detid);
        isHalfCell[i] = rhtools_.isHalfCell(detid);
      }
    }

    auto recHitsTable =
        std::make_unique<nanoaod::FlatTable>(nRecHits, tableName_, /*singleton*/ false, /*extension*/ false);
    recHitsTable->addColumn<float>("x", posX, "RecHits position x [cm]", precision_);
    recHitsTable->addColumn<float>("y", posY, "RecHits position Y [cm]", precision_);
    recHitsTable->addColumn<float>("Z", posZ, "RecHits position Z [cm]", precision_);
    recHitsTable->addColumn<int>("layer", layer, "RecHits layer number", precision_);
    recHitsTable->addColumn<float>("energy", energy, "RecHits energy [GeV]", precision_);
    recHitsTable->addColumn<float>("thickness", thickness, "RecHits wafer thickness, 0 means scintillator", precision_);
    recHitsTable->addColumn<bool>("isHalfCell", isHalfCell, "RecHits is half cell", precision_);
    iEvent.put(std::move(recHitsTable), tableName_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<bool>("skipNonExistingSrc", false)
        ->setComment("whether or not to skip producing the table on absent input product");
    desc.add<std::string>("tableName", "hltMergeRecHits")->setComment("name of the flat table ouput");
    desc.add<std::vector<edm::InputTag>>("recHits", {edm::InputTag("HGCRecHit", "HGCEERecHits")});
    desc.add<int>("precision", 7);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const bool skipNonExistingSrc_;
  const std::string tableName_;
  const unsigned int precision_;
  std::vector<edm::EDGetTokenT<std::vector<HGCRecHit>>> hits_token_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geom_token_;
  hgcal::RecHitTools rhtools_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecHitsExtraTableProducer);
