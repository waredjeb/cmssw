#include <algorithm>
#include <set>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TrackstersLinkingbyPCA.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"

using namespace ticl;
TrackstersLinkingbyPCA::TrackstersLinkingbyPCA(const edm::ParameterSet &conf, const CacheBase *cache)
    : TrackstersLinkingAlgoBaseT(conf, cache){};

TrackstersLinkingbyPCA::~TrackstersLinkingbyPCA(){};
void TrackstersLinkingbyPCA::trackstersInfo(const typename TrackstersLinkingAlgoBaseT::Inputs &input) {
  const auto tracksters = input.tracksters;
  for (auto t : tracksters) {
    auto e_over_h = (t.raw_em_pt() / ((t.raw_pt() - t.raw_em_pt()) != 0. ? (t.raw_pt() - t.raw_em_pt()) : 1.));
    std::cout << "\nTrackster raw_pt: " << t.raw_pt() << " raw_em_pt: " << t.raw_em_pt() << " eoh: " << e_over_h
              << " barycenter: " << t.barycenter() << " eta,phi (baricenter): " << t.barycenter().eta() << ", "
              << t.barycenter().phi() << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", "
              << t.eigenvectors(0).phi()
              << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy()
              << " seedID: " << t.seedID() << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size()
              << " average usage: "
              << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
                  (float)t.vertex_multiplicity().size())
              << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
              << " probs(ga/e/mu/np/cp/nh/am/unk): ";
    for (auto const &p : t.id_probabilities()) {
      std::cout << std::fixed << p << " ";
    }
    std::cout << " sigmas: ";
    for (auto const &s : t.sigmas()) {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
}
