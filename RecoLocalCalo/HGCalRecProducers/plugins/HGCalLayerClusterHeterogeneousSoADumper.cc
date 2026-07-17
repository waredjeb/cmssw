#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <fmt/format.h>

class HGCalLayerClusterHeterogeneousSoADumper : public edm::global::EDAnalyzer<> {
public:
  HGCalLayerClusterHeterogeneousSoADumper(edm::ParameterSet const& iConfig)
      : token_{consumes(iConfig.getParameter<edm::InputTag>("src"))} {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hltHgcalSoALayerClustersProducer"));
    descriptions.addWithDefaultLabel(desc);
  }

  void analyze(edm::StreamID iStream, edm::Event const& iEvent, edm::EventSetup const& iSetup) const override {
    auto const& data = iEvent.get(token_);

    auto const view = data.view();
    auto const position_v = view.position();
    auto const energy_v = view.energy();
    const int numberOfClusters = position_v.metadata().size();
    std::cout << fmt::format("hgcalSoALayerClustersProducer size = {}", numberOfClusters) << std::endl;
    for (int i = 0; i < numberOfClusters; ++i) {
      std::cout << fmt::format("CLUSTERS_SOA {}, energy = {:.{}f}, x = {:.{}f}, y = {:.{}f}, z= {:.{}f}",
                               i,
                               energy_v[i].energy(),
                               std::numeric_limits<float>::max_digits10,
                               position_v[i].x(),
                               std::numeric_limits<float>::max_digits10,
                               position_v[i].y(),
                               std::numeric_limits<float>::max_digits10,
                               position_v[i].z(),
                               std::numeric_limits<float>::max_digits10)
                << std::endl;
    }
  }

private:
  edm::EDGetTokenT<reco::CaloClusterHostCollection> const token_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClusterHeterogeneousSoADumper);
