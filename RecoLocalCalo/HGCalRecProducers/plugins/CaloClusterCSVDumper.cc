// Minimal analyzer to dump reco::CaloCluster collections to CSV, for comparing
// the layer-clustering engines (legacy CLUE / alpaka CLUEstering / old device CLUE)
// across offline and HLT. One row per cluster.
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include <fstream>
#include <mutex>

class CaloClusterCSVDumper : public edm::one::EDAnalyzer<> {
public:
  explicit CaloClusterCSVDumper(edm::ParameterSet const& ps)
      : token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("src"))),
        tag_(ps.getParameter<std::string>("tag")),
        fileName_(ps.getParameter<std::string>("out")) {
    std::ofstream f(fileName_, std::ios::out | std::ios::trunc);
    f << "tag,run,lumi,event,icl,seeddet,seedid,energy,eta,phi,x,y,z,nhits\n";
  }

  void analyze(edm::Event const& e, edm::EventSetup const&) override {
    auto const& clusters = e.get(token_);
    std::lock_guard<std::mutex> lock(mtx_);
    std::ofstream f(fileName_, std::ios::out | std::ios::app);
    int icl = 0;
    for (auto const& c : clusters) {
      f << tag_ << ',' << e.id().run() << ',' << e.id().luminosityBlock() << ',' << e.id().event() << ',' << icl++
        << ',' << c.seed().det()  // Detector enum: HGCalEE=8, HGCalHSi=9, HGCalHSc=10 -> silicon = {8,9}
        << ',' << c.seed().rawId()  // seed DetId raw id -> key for seed-by-seed matching
        << ',' << c.energy() << ',' << c.eta() << ',' << c.phi() << ',' << c.x() << ',' << c.y() << ',' << c.z()
        << ',' << c.hitsAndFractions().size()
        << '\n';
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& d) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hgcalMergeLayerClusters"));
    desc.add<std::string>("tag", "algo");
    desc.add<std::string>("out", "clusters.csv");
    d.addWithDefaultLabel(desc);
  }

private:
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> token_;
  std::string tag_;
  std::string fileName_;
  std::mutex mtx_;
};

DEFINE_FWK_MODULE(CaloClusterCSVDumper);
