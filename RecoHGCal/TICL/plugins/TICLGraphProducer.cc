#include <memory>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

using namespace ticl;
typedef math::XYZVector Vector;
class TICLGraphProducer : public edm::stream::EDProducer<> {
public:
  explicit TICLGraphProducer(const edm::ParameterSet &ps);
  ~TICLGraphProducer() override{};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  void initialize(const HGCalDDDConstants *hgcons,
                  const hgcal::RecHitTools rhtools,
                  const edm::ESHandle<MagneticField> bfieldH,
                  const edm::ESHandle<Propagator> propH);
  void buildLayers();
  void beginJob();
  void endJob();

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;

private:
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const StringCutObjectSelector<reco::Track> cutTk_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  const std::string detector_;
  const std::string propName_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  const double track_sep_;
  const double trackster_sep_;
  const double delta_etaphi_;

  const HGCalDDDConstants *hgcons_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  edm::ESHandle<MagneticField> bfield_;
  edm::ESHandle<Propagator> propagator_;
  std::unique_ptr<GeomDet> firstDisk_[2];
};

TICLGraphProducer::TICLGraphProducer(const edm::ParameterSet &ps)
    : tracksters_clue3d_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      tracks_token_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("tracks"))),
      cutTk_(ps.getParameter<std::string>("cutTk")),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      detector_(ps.getParameter<std::string>("detector")),
      propName_(ps.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(
          esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(edm::ESInputTag("", propName_))),
      track_sep_(ps.getParameter<double>("trackSep")),
      trackster_sep_(ps.getParameter<double>("tracksterSep")),
      delta_etaphi_(ps.getParameter<double>("deltaEtaPhi")) {
  produces<std::vector<TICLGraph>>();
//  produces<std::vector<TICLGraph>>("fromTracksters");
  std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ =
      esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", detectorName_));
}

void TICLGraphProducer::beginJob() {}

void TICLGraphProducer::endJob(){};

void TICLGraphProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();

  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);

  bfield_ = es.getHandle(bfield_token_);
  propagator_ = es.getHandle(propagator_token_);
  buildLayers();
};

void TICLGraphProducer::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;

    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());
  }
}

float separation(const Trackster &t1, const Trackster &t2) {
  auto const &bary1 = t1.barycenter();
  auto const &bary2 = t2.barycenter();

  auto r1 = std::sqrt(bary1.x() * bary1.x() + bary1.y() * bary1.y());
  auto const r1OverZ1 = r1 / std::abs(bary1.z());
  auto r2 = std::sqrt(bary2.x() * bary2.x() + bary2.y() * bary2.y());
  auto r2overz2 = r2 / std::abs(bary2.z());
  r2 = r2overz2 * bary1.z();
  r1 = r1OverZ1 * bary1.z();
  auto delta_phi = reco::deltaPhi(bary1.phi(), bary2.phi());
  //std::cout << "r1 " << r1 << " r2 " << r2 << " z1 " << bary1.z() << " deltaPhi " << delta_phi << std::endl;
  return std::sqrt((r1 - r2) * (r1 - r2) + r2 * r2 * delta_phi * delta_phi);
}
float updateNode(Node &node,
                 std::vector<Node> &nodes,
                 std::vector<Trackster> &tracksters,
                 uint32_t t2,
                 const float sep_th,
                 std::vector<int> &mask_tracksters_for_node,
                 std::vector<int> &mask_tracksters_for_edge,
                 std::unordered_map<uint32_t, uint32_t> trackster_to_node,
                 std::string &tabs) {
  auto const &t1 = node.getId();
  auto const &trackster1 = tracksters[t1];
  auto const &trackster2 = tracksters[t2];
  auto distance = separation(trackster1, trackster2);
  auto dotProduct = [](const Vector &v1, const Vector &v2) {
    return std::abs(v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z());
  };
  float score = -1.f;
  mask_tracksters_for_node[t1] = 1;
  //  std::cout << tabs << "Masking " << t1 << std::endl;
  //  std::cout << tabs << "Update node " << t1 << " with " << t2 << " distance " << distance  << " th " << sep_th << std::endl;
  if (distance < sep_th && mask_tracksters_for_edge[t2] == 0) {
//    std::cout << "Trackster 1 " << trackster1.barycenter() << " Trackster 2 " << trackster2.barycenter() << std::endl;
    auto const &dir1 = trackster1.eigenvectors(0);
    auto const &dir2 = trackster2.eigenvectors(0);
    auto dot = dotProduct(dir1, dir2);
    score = 1 - (distance / sep_th + (1 - dot)) / 2;
    //  std::cout << "Adding edge between " << t1 << " and " << t2 << " with score " << score << std::endl;
    node.addEdge(t2, score);
    auto otherNode = std::find_if(nodes.begin(), nodes.end(), [=](Node &n) { return n.getId() == t2 ;});
    bool found = false;
    if (otherNode != nodes.end()) {
      for (auto const &[neigh, weigh] : otherNode->getWeightedEdges()) {
        if (neigh == node.getId()) {
          found = true;
        }
      }
      if (!found) {
        otherNode->addEdge(node.getId(), score);
      }
    }
  }
  return score;
}

void buildGraph(std::vector<Node> &nodes,
                Node &node,
                std::array<TICLLayerTile, 2> &trackster_tiles,
                std::vector<Trackster> &tracksters,
                const uint32_t tracksterId,
                const Node &lastNode,
                const float sep_th,
                const float delta,
                std::vector<int> &mask_tracksters_for_node,
                std::vector<int> &mask_tracksters_for_edge,
                std::unordered_map<uint32_t, uint32_t> trackster_to_node,
                std::string &tabs) {
  assert(node.getId() == tracksterId);
  if (mask_tracksters_for_node[tracksterId] == 0) {
    //    std::cout << "Starting building graph for node " << node.getId() << std::endl;
    auto const &trackster = tracksters[tracksterId];
    auto const &barycenter = trackster.barycenter();
    int sideZ = barycenter.eta() > 0;
    auto const &tile = trackster_tiles[sideZ];
    auto tracksterRoverZ =
        std::sqrt(barycenter.x() * barycenter.x() + barycenter.y() * barycenter.y()) / std::abs(barycenter.z());

    double eta_min = std::max(abs(barycenter.eta()) - delta, (double)TileConstants::minEta);
    double eta_max = std::min(abs(barycenter.eta()) + delta, (double)TileConstants::maxEta);

    std::array<int, 4> search_box =
        tile.searchBoxEtaPhi(eta_min, eta_max, barycenter.phi() - delta, barycenter.phi() + delta);
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
        for (const unsigned &t_i : in_tile) {
          //          std::cout << tabs << "Node " << tracksterId <<  " in Tile " << t_i <<  " LastNode " << lastNode.getId() << std::endl;
          if (t_i != tracksterId &&
              ((t_i != lastNode.getId() && lastNode.isTrackster() == 1) || lastNode.isTrackster() != 1)) {
            //update graph!
            auto score = updateNode(node, nodes, tracksters, t_i, sep_th, mask_tracksters_for_node, mask_tracksters_for_edge, trackster_to_node, tabs);
            //           std::cout << tabs << "Updating node " << score << std::endl;
            if (score >= 0) {
              //           std::cout << tabs << "Creating new node " << t_i << std::endl;
              Node newNode(t_i, 1);

              //              std::cout << tabs << "Adding edge between " <<  t_i << " and " << tracksterId << " with score " << score << std::endl;
              newNode.addEdge(tracksterId, score);
              tabs += "\t";
              buildGraph(nodes,
                         newNode,
                         trackster_tiles,
                         tracksters,
                         t_i,
                         node,
                         sep_th,
                         delta,
                         mask_tracksters_for_node,
                         mask_tracksters_for_edge,
                         trackster_to_node,
                         tabs);
            }
          }
        }
      }
    }
    tabs.pop_back();
    nodes.push_back(node);
  }
}

void TICLGraphProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  //  std::cout << "TICL GRAPH! " << std::endl;
  edm::Handle<std::vector<Trackster>> trackstersclue3d_h;
  evt.getByToken(tracksters_clue3d_token_, trackstersclue3d_h);
  auto trackstersclue3d = *trackstersclue3d_h;

  auto const &tracks = evt.get(tracks_token_);

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);
  //std::vector<Trackster> trackstersclue3d_sorted(trackstersclue3d);
  //std::sort(trackstersclue3d_sorted.begin(), trackstersclue3d_sorted.end(), [](Trackster& t1, Trackster& t2){return t1.barycenter().z() < t2.barycenter().z();});

  //Fill tiles
  TICLLayerTile tracksterTilePos;
  TICLLayerTile tracksterTileNeg;
  std::vector<int> mask_tracksters_for_node(trackstersclue3d.size(), 0);
  std::vector<int> mask_tracksters_for_edge(trackstersclue3d.size(), 0);

  auto distance = [](float r0, float r1, float phi0, float phi1) {
    auto delta_phi = reco::deltaPhi(phi0, phi1);
    return std::sqrt((r0 - r1) * (r0 - r1) + r1 * r1 * delta_phi * delta_phi);
  };

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];
    if (t.barycenter().eta() > 0.) {
      tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    } else if (t.barycenter().eta() < 0.) {
      tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    }
  }

  std::array<TICLLayerTile, 2> tiles = {{tracksterTileNeg, tracksterTilePos}};
  std::vector<Node> allNodes;
  std::vector<Node> trackNodes;
  trackNodes.reserve(tracks.size());
  auto graphs = std::make_unique<std::vector<TICLGraph>>();
  for (size_t i_track = 0; i_track < tracks.size(); ++i_track) {
    std::vector<std::pair<uint32_t, float>> scores;
    auto const &track = tracks[i_track];
    if (cutTk_(track)) {
      Node trackNode(i_track, 0);
      int iSide = int(track.outerEta() > 0);
//      std::cout << " iSide " << iSide << std::endl;
      auto tile = tiles[iSide];
      const auto &fts = trajectoryStateTransform::outerFreeState((track), bFieldProd);
      const auto &tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      const auto &tsosPosition = tsos.globalPosition();
      const auto &tsosDirection = tsos.globalMomentum().unit();
      const auto tsosZ = tsosPosition.z();
      const auto tsosPhi = tsosPosition.phi();
      auto trackRoverZ =
          std::sqrt(tsosPosition.x() * tsosPosition.x() + tsosPosition.y() * tsosPosition.y()) / std::abs(tsosZ);

      if (tsos.isValid()) {
        auto const trackPhi = tsosPosition.phi();
        auto const trackEta = tsosPosition.eta();
        double eta_min = std::max(abs(trackEta) - delta_etaphi_, (double)TileConstants::minEta);
        double eta_max = std::min(abs(trackEta) + delta_etaphi_, (double)TileConstants::maxEta);

        std::array<int, 4> search_box =
            tile.searchBoxEtaPhi(eta_min, eta_max, trackPhi - delta_etaphi_, trackPhi + delta_etaphi_);

        for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
          for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
            const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];

            for (const unsigned &t_i : in_tile) {
              if (mask_tracksters_for_node[t_i] == 0) {
                auto const &trackster = trackstersclue3d[t_i];
                auto const &tracksterBarycenter = trackster.barycenter();
                auto const tracksterPhi = tracksterBarycenter.phi();
                auto tracksterRoverZ = std::sqrt(tracksterBarycenter.x() * tracksterBarycenter.x() +
                                                 tracksterBarycenter.y() * tracksterBarycenter.y()) /
                                       std::abs(tracksterBarycenter.z());
                const auto sep = distance(trackRoverZ * tsosZ, tracksterRoverZ * tsosZ, tsosPhi, tracksterPhi);
                if (sep <= track_sep_) {
                  auto tracksterEigenVector = trackster.eigenvectors(0);
                  auto tracksterDir =
                      Global3DVector(tracksterEigenVector.x(), tracksterEigenVector.y(), tracksterEigenVector.z());
                  auto const dotProduct = tsosDirection.dot(tracksterDir);
                  auto const score = 1 - (sep / track_sep_ + (1 - dotProduct)) / 2;
                  scores.emplace_back(t_i, score);
                }
              }
            }
          }
        }
        //get best score trackster
        if (!scores.empty()) {
          auto max = std::max_element(
              scores.begin(), scores.end(), [](std::pair<uint32_t, float> v1, std::pair<uint32_t, float> v2) {
                return v1.second < v2.second;
              });
          auto argmaxVal = std::distance(scores.begin(), max);
          auto const bestTracksterId = scores[argmaxVal].first;
          auto const bestScore = scores[argmaxVal].second;
          trackNodes.push_back(trackNode);
          Node tracksterNode(bestTracksterId, 1);
          trackNode.addEdge(bestTracksterId, bestScore);
          std::vector<Node> graphFromTrack;
          std::string tabs = "";
          std::unordered_map<uint32_t, uint32_t> trackster_to_node;
          buildGraph(graphFromTrack,
                     tracksterNode,
                     tiles,
                     trackstersclue3d,
                     tracksterNode.getId(),
                     trackNode,
                     trackster_sep_,
                     delta_etaphi_,
                     mask_tracksters_for_node,
                     mask_tracksters_for_edge,
                     trackster_to_node,
                     tabs);
          mask_tracksters_for_node[bestTracksterId] = 1;
          mask_tracksters_for_edge = mask_tracksters_for_node;
          auto g = TICLGraph(graphFromTrack);
          g.setTrackToTracksterEdge(i_track, bestTracksterId);
          graphs->push_back(g);
        }
      }
    }
  }

//  for (size_t ig = 0; ig < graphs->size(); ++ig) {
//    std::cout << "@@@@ Track Graph @@@@ " << ig << std::endl;
//    for (size_t in = 0; in < (*graphs)[ig].size(); ++in) {
//      auto nod = (*graphs)[ig].getNode(in);
//        std::cout << "\t Node " << in << " Trackster " << nod.getId() << " Mask " << mask_tracksters_for_node[nod.getId()] << std::endl;
//      auto edges = nod.getWeightedEdges();
//      for (auto const &edge : edges) {
//           std::cout << "\t\t" << " Trackster " << edge.first << " Score " << edge.second << std::endl;
//      }
//    }
//  }

  auto tracksterGraphs = std::make_unique<std::vector<TICLGraph>>();
 // std::cout << "#################### BUILDING TRACKSTER NODES ####################" << std::endl;
  for (size_t i = 0; i < mask_tracksters_for_node.size(); ++i) {
    if (mask_tracksters_for_node[i] == 0) {
      std::vector<Node> tracksterNodes;
      Node tracksterNode(i, 1);
      std::string tabs = "";
      std::unordered_map<uint32_t, uint32_t> trackster_to_node;

      buildGraph(tracksterNodes,
                 tracksterNode,
                 tiles,
                 trackstersclue3d,
                 tracksterNode.getId(),
                 tracksterNode,
                 trackster_sep_,
                 delta_etaphi_,
                 mask_tracksters_for_node,
                 mask_tracksters_for_edge,
                 trackster_to_node,
                 tabs);
      mask_tracksters_for_edge = mask_tracksters_for_node;
      tracksterGraphs->emplace_back(tracksterNodes);
    }
  }
//  for (size_t ig = 0; ig < tracksterGraphs->size(); ++ig) {
//       std::cout << "#### Tracksters Graph #### " << ig << std::endl;
//    for (size_t in = 0; in < (*tracksterGraphs)[ig].size(); ++in) {
//      auto nod = (*tracksterGraphs)[ig].getNode(in);
//        std::cout << "\t Node " << in << " Trackster " << nod.getId() << std::endl;
//      auto edges = nod.getWeightedEdges();
//      for (auto const &edge : edges) {
//             std::cout << "\t\t" << " Trackster " << edge.first << " Score " << edge.second << std::endl;
//      }
//    }
//  }
  //  auto resultGraph = std::make_unique<TICLGraph>(allNodes);
  //
  
  std::copy(std::begin(*tracksterGraphs), std::end(*tracksterGraphs), std::back_inserter(*graphs));
  evt.put(std::move(graphs));
  //evt.put(std::move(graphs), "fromTracks");
  //evt.put(std::move(tracksterGraphs), "fromTracksters");
}

void TICLGraphProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons1stStep"));
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("trackSep", 10.f);
  desc.add<double>("tracksterSep", 10.f);
  desc.add<double>("deltaEtaPhi", 0.05);
  descriptions.add("ticlGraphProducer", desc);
}

DEFINE_FWK_MODULE(TICLGraphProducer);
