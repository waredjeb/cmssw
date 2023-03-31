#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByLouvainAlgo.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TrackstersPCA.h"
using namespace ticl;

LinkingAlgoByLouvainAlgo::LinkingAlgoByLouvainAlgo(const edm::ParameterSet &conf)
    : LinkingAlgoBase(conf),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      del_ts_em_had_(conf.getParameter<double>("delta_ts_em_had")),
      del_ts_had_had_(conf.getParameter<double>("delta_ts_had_had")),
      separationSmall_threshold_(conf.getParameter<double>("separationSmall")),
      separation_threshold_(conf.getParameter<double>("separation")),
      maxDepth_(conf.getParameter<int>("maxDepth")),
      timing_quality_threshold_(conf.getParameter<double>("track_time_quality_threshold")),
      cutTk_(conf.getParameter<std::string>("cutTk")) {}

LinkingAlgoByLouvainAlgo::~LinkingAlgoByLouvainAlgo() {}

void LinkingAlgoByLouvainAlgo::initialize(const HGCalDDDConstants *hgcons,
                                          const hgcal::RecHitTools rhtools,
                                          const edm::ESHandle<MagneticField> bfieldH,
                                          const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}

void LinkingAlgoByLouvainAlgo::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide),
                    Disk::RotationType(),
                    SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
            .get());
  }
}
void LinkingAlgoByLouvainAlgo::dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection,
                                              const char *label) const {
  //#ifdef EDM_ML_DEBUG
  //  if (!(LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced))
  //    return;

  ////std::cout << "All links found - " << label << "\n";
  ////std::cout << "(seed can either be a track or trackster depending on the step)\n";
  for (unsigned i = 0; i < resultCollection.size(); ++i) {
    ////std::cout << "seed " << i << " - tracksters : ";
    const auto &links = resultCollection[i];
    for (unsigned j = 0; j < links.size(); ++j) {
      ////std::cout << j;
    }
    ////std::cout << "\n";
  }
  //#endif  // EDM_ML_DEBUG
}

class Community {
public:
  Community(const TICLGraph &g, const std::vector<uint32_t> &nodes, float iw, float tw)
      : g_(g), nodes_(nodes), internal_weight_(iw), total_weight_(tw){};
  void updateTotalWeight() {
    float community_total_weight = 0;
    for (size_t j = 0; j < nodes_.size(); j++) {
      const Node &node_data = g_.getNode(nodes_[j]);
      community_total_weight += node_data.getWeightedDegree();
    }
    total_weight_ = community_total_weight;
  }
  bool isNodeInCommunity(uint32_t node_id) const {
    return std::find(nodes_.begin(), nodes_.end(), node_id) != nodes_.end();
  }
  void updateInternalWeight() {
    float tmp_internal_weight = 0;
    for (size_t i = 0; i < nodes_.size(); i++) {
      const Node &node_data = g_.getNode(nodes_[i]);
      for (const auto &[neighbor, weight] : node_data.getWeightedEdges()) {
        if (isNodeInCommunity(neighbor)) {
          tmp_internal_weight += weight;
        }
      }
    }
    internal_weight_ = tmp_internal_weight;
  }

  std::vector<uint32_t> getNodes() const { return nodes_; }
  float getInternalWeight() const { return internal_weight_; }
  float getTotalWeight() const { return total_weight_; }
  const TICLGraph &getGraph() const { return g_; }
  void setInternalWeight(const float w) { internal_weight_ = w; }
  void setTotalWeight(const float w) { total_weight_ = w; }
  void addInternalWeight(const float w) { internal_weight_ += w; }
  void addTotalWeight(const float w) { total_weight_ += w; }

private:
  const TICLGraph &g_;
  std::vector<uint32_t> nodes_;
  float internal_weight_;
  float total_weight_;
};
std::unordered_map<uint32_t, uint32_t> map_trackster_to_node(const TICLGraph &graph) {
  std::unordered_map<uint32_t, uint32_t> result;
  for (size_t i = 0; i < graph.size(); i++) {
    auto node = graph.getNode(i);
    if (node.isTrackster())
      result[node.getId()] = i;
  }
  return result;
}
std::unordered_map<uint32_t, uint32_t> map_node_to_community(const std::vector<Community> &communities) {
  std::unordered_map<uint32_t, uint32_t> result;
  for (size_t i = 0; i != communities.size(); i++) {
    for (int node : communities[i].getNodes()) {
      result[node] = i;
    }
  }
  return result;
}
float modularity(const TICLGraph &graph, std::vector<Community> &communities) {
  float q = 0.f;
  const float m2 = graph.getTotalWeight();
  for (size_t i = 0; i < communities.size(); ++i) {
    auto const &comm = communities[i];
    q += comm.getInternalWeight() / m2 - ((comm.getTotalWeight() / m2) * (comm.getTotalWeight() / m2));
  }
  return q;
}

float getWeightedEdgesInCommunity(const Node &node,
                                  const int communityIndex,
                                  std::unordered_map<uint32_t, uint32_t> &node_to_community,
                                  std::unordered_map<uint32_t, uint32_t> &trackster_to_node) {
  //    std::unordered_map<uint32_t, float> weighted_edges;
  float result = 0.f;
  for (const auto &[neighbor, weight] : node.getWeightedEdges()) {
    auto neighbor_i = trackster_to_node[neighbor];
    int neighbor_community = node_to_community.at(neighbor_i);
    if (neighbor_community == communityIndex) {
      result += weight;
    }
  }
  return result;
}
void removeNodeFromComm(const Node &node,
                        const uint32_t node_i,
                        Community &community,
                        const float weightInCommunity,
                        std::unordered_map<uint32_t, uint32_t> &n2c) {
  community.addTotalWeight(-node.getWeightedDegree());
  //  community.total_weight -= node.getWeightedDegree();
  community.addInternalWeight(-2 * weightInCommunity + community.getGraph().getSelfWeight(node_i));
  //  community.internal_weight -= 2 * weightInCommunity + graph.getSelfWeight(node_i);
  n2c[node_i] = -1;
}

float modularityGain(
    const Node &node, const uint32_t node_i, Community &comm, const uint32_t comm_i, const float weightInCommunity) {
  float totc = comm.getTotalWeight();
  float degc = comm.getInternalWeight();
  float m2 = comm.getGraph().getTotalWeight();

  return (weightInCommunity - totc * degc / m2);
}

void insertNodeInComm(const Node &node,
                      const uint32_t node_i,
                      Community &community,
                      const uint32_t comm_i,
                      const float weightInCommunity,
                      std::unordered_map<uint32_t, uint32_t> &n2c) {
  community.addTotalWeight(node.getWeightedDegree());
  //  community.total_weight += node.getWeightedDegree();
  community.addInternalWeight(2 * weightInCommunity + community.getGraph().getSelfWeight(node_i));
  //  community.internal_weight += 2 * weightInCommunity + graph.getSelfWeight(node_i);
  n2c[node_i] = comm_i;
}

bool one_level(const TICLGraph &graph,
               std::vector<Community> &communities,
               std::unordered_map<uint32_t, uint32_t> &trackster_to_node,
               std::unordered_map<uint32_t, uint32_t> &node_to_community) {
  /*
  graph.getNode(i) returns node i-th in the graph
  node.getWeightedDegree() returns the total degree of the node
  community.updateTotalWeight(graph) updates the total degree of the community
  community.updateInternalWeight(graph) updates the internal degree of the community
  map_node_to_community(communities) returns std::unordered_map between a node and its community
  trackster_to_node is a map between the trackster Id and the node id (id of the node in the graph)
  getWeightedEdgesInCommunity(node,i,node_to_community, trackster_to_node) returns only the contribution of node to the internal degree of the community
  graph.getTotalWeight() returns the total degree of the graph (already double counting the edges since the graph is undirected)
  */
  assert(communities.size() == graph.size());
  std::unordered_map<uint32_t, uint32_t> new_node_to_community;
  bool improvement = false;
  auto modImprov = true;
  auto new_mod = modularity(graph, communities);
  auto cur_mod = new_mod;
  auto node_move = 1;
  auto iterations = 0;
  while (node_move >= 1 && modImprov) {
    node_move = 0;
    modImprov = false;

    for (size_t node_i = 0; node_i < graph.size(); node_i++) {
      auto comm = node_to_community[node_i];
      auto &node = graph.getNode(node_i);
      auto node_w = node.getWeightedDegree();
      auto degInComm = getWeightedEdgesInCommunity(node, comm, node_to_community, trackster_to_node);
      //      std::cout << "Node " << node_i << " Weight " << node_w << " Weight in Comm " << comm << " " << degInComm
      //              << std::endl;
      removeNodeFromComm(node, node_i, communities[comm], degInComm, node_to_community);
      auto best_comm = comm;
      auto best_nblinks = 0.f;
      auto best_increase = 0.f;
      cur_mod = new_mod;
      for (auto &[neigh, n_w] : node.getWeightedEdges()) {
        auto comm_n = node_to_community[trackster_to_node[neigh]];
        auto dnc = getWeightedEdgesInCommunity(node, comm_n, node_to_community, trackster_to_node);
        //std::cout << "trackster to node " << trackster_to_node[neigh] << " trackster " << neigh << std::endl;
        //std::cout << "node " << node_i << " comm n " << comm_n << std::endl;
        auto commN = communities[comm_n];
        auto increase = modularityGain(node, node_i, communities[comm_n], comm_n, dnc);
        //        std::cout << "\t Neighbour " << trackster_to_node[neigh] << " community " << comm_n << " gain mod " << increase
        //                << std::endl;
        if (increase > best_increase) {
          best_comm = comm_n;
          best_nblinks = getWeightedEdgesInCommunity(node, comm_n, node_to_community, trackster_to_node);
          //         std::cout << " Weight in comm " << comm_n << " " << best_nblinks << std::endl;
          best_increase = increase;
        }
      }
      //   std::cout << "Best comm " << best_comm << " N2C " << comm << std::endl;
      insertNodeInComm(node, node_i, communities[best_comm], best_comm, best_nblinks, node_to_community);
      //     std::cout << "After Best comm " << best_comm << " N2C " << node_to_community[node_i] << std::endl;
      if (best_comm != comm)
        node_move++;

      float total_tot = 0.f;
      float total_in = 0.f;
      for (size_t i_comm = 0; i_comm < communities.size(); i_comm++) {
        auto const &c = communities[i_comm];
        total_tot += c.getTotalWeight();
        total_in += c.getInternalWeight();
      }
      new_mod = modularity(graph, communities);
      //     std::cout << "New mod " << new_mod << " Cur Mod " << cur_mod << " diff " << std::abs(new_mod - cur_mod)
      //             << " Node moved " << node_move << std::endl;
      if (std::abs(new_mod - cur_mod) > 10e-5) {
        modImprov = true;
      }
    }
    if (node_move >= 1) {
      improvement = true;
    }
    auto s = improvement == true ? "True" : "False";
    //   std::cout << "Iterations " << iterations << " Node move " << node_move << " improvement " << s << std::endl;
    iterations++;
  }
  return improvement;
}

std::vector<Community> louvain(const TICLGraph &graph,
                               std::unordered_map<uint32_t, uint32_t> &trackster_to_node,
                               std::vector<Community> &communities,
                               std::unordered_map<uint32_t, uint32_t> &node_to_community) {
  auto improvement = true;
//  auto mod = modularity(graph, communities);
  //float new_mod;
 // int level = 0;
  while (improvement) {
    improvement = false;
    improvement = one_level(graph, communities, trackster_to_node, node_to_community);
//    new_mod = modularity(graph, communities);
//    std::cout << "Modularity increased from " << mod << " To " << new_mod << std::endl;
 //   mod = new_mod;
//    level++;
  }
//  std::cout << "Level " << level << std::endl;
  return communities;
}

void LinkingAlgoByLouvainAlgo::linkTracksters(const std::vector<TICLGraph> &graphsFromTrack,
                                           //   const std::vector<TICLGraph> &graphFromTracksters,
                                              const edm::Handle<std::vector<reco::Track>> tkH,
                                              const edm::ValueMap<float> &tkTime,
                                              const edm::ValueMap<float> &tkTimeErr,
                                              const edm::ValueMap<float> &tkTimeQual,
                                              const std::vector<reco::Muon> &muons,
                                              const edm::Handle<std::vector<Trackster>> tsH,
                                              const std::vector<reco::CaloCluster> &layerClusters,
                                              const edm::ValueMap<std::pair<float, float>> &layerClustersTimes,
                                              std::vector<Trackster> &resultTrackstersMerged,
                                              std::vector<TICLCandidate> &candidates,
                                              std::vector<TICLCandidate> &chargedCandidatesFromTracks,
                                              const EnergyRegressionAndIDModel &model) {
  //std::cout << "LOUVAIN ALGO " << std::endl;
  const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);
  //needed to keep track of trackster indices when building TICLCandidate
  std::vector<std::vector<unsigned>> tracksterMergeSmallCollectionIndices;
  std::vector<std::vector<unsigned>> tracksterMergeCollectionIndices;
  //needed to keep track of trackster indices when building TICLCandidate
  std::vector<Trackster> tracksterMergeSmallCollection;
  std::vector<Trackster> tracksterMergeCollection;
  //after cleaning needed to keep track of trackster indices when building TICLCandidate
  std::vector<Trackster> resultTrackstersSmallMerged;  //after cleaning.
  std::vector<std::vector<unsigned>> resultTrackstersSmallMergedIndices;
  std::vector<std::vector<unsigned>> resultTrackstersMergedIndices;

  std::array<TICLLayerTile, 2> tracksterPropTiles = {};       // all Tracksters
  std::array<TICLLayerTile, 2> tracksterSmallPropTiles = {};  // all Tracksters
  std::array<TICLLayerTile, 2> tracksPropTiles = {};          // all Tracks

  if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
    LogDebug("LinkingAlgoByLouvainAlgo") << "------- Geometric Linking ------- \n";

  //std::cout << "Graphs from tracks " << graphsFromTrack.size() << std::endl;
  //std::cout << "Starting Louvain method on all the graphs from tracks " << std::endl;

  std::vector<std::vector<std::vector<uint32_t>>> g_communities;
  auto g_i = 0;

  auto graphs = graphsFromTrack;
//  std::copy(std::begin(graphFromTracksters), std::end(graphFromTracksters), std::back_inserter(graphs));

  std::vector<std::vector<bool>> chargedCommunities;
  for (auto const &g : graphs) {
    //initialize communities, one community for each node
    std::vector<Community> communities;
    for (size_t i = 0; i != g.size(); i++) {
      const float initial_weight = g.getNode(i).getWeightedDegree();
      std::vector<uint32_t> n = {{static_cast<uint32_t>(i)}};
      communities.emplace_back(g, n, 0., 0.);
    }
    //create map between node and community
    auto node_to_community = map_node_to_community(communities);
    //create map between trackster id and node id
    auto trackster_to_node = map_trackster_to_node(g);
    //start louvain algorithm on graph g
    if (g.size() > 1) {
      louvain(g, trackster_to_node, communities, node_to_community);
    }
    std::vector<std::vector<uint32_t>> r_communities;
    std::vector<bool> chargedCommunitiesInGraph;

    //collect nodes for each community
    // Iterate over all the nodes and assign each node to its corresponding community
    for (auto const &[node, community] : node_to_community) {
      // Check if the community vector exists, otherwise create it
      if (community >= r_communities.size()) {
        r_communities.resize(community + 1);
        chargedCommunitiesInGraph.resize(community + 1);
      }
      // Add the node to its corresponding community
      r_communities[community].push_back(node);
      if (static_cast<int>(node) == g.getTrackToTracksterEdge().second) {
        chargedCommunitiesInGraph.push_back(true);
      } else {
        chargedCommunitiesInGraph.push_back(false);
      }
    }
    chargedCommunities.push_back(chargedCommunitiesInGraph);
    g_communities.push_back(r_communities);
    g_i++;
  }

  //build trackster merged and ticl candidates
  auto i_c = 0;
  for (auto const &g_c : g_communities) {
//    std::cout << "Printing community for graph " << i_c << std::endl;
    auto i_cc = 0;
  //  std::cout << "Community size " <<  g_c.size() << std::endl;
     for (auto const &c : g_c) {
    //  std::cout << "Community " << i_cc << std::endl;
      if (c.size() > 0) {
        Trackster outTrackster;
        auto updatedSize = outTrackster.vertices().size();
        TICLCandidate candidate;
        for (auto const &n : c) {
    //      std::cout << "\t"
      //              << " Node " << n << std::endl;
          auto trackster_id = graphs[i_c].getNode(n).getId();

          auto const &thisTrackster = tracksters[trackster_id];
          updatedSize += thisTrackster.vertices().size();
          outTrackster.vertices().reserve(updatedSize);
          outTrackster.vertex_multiplicity().reserve(updatedSize);
          std::copy(std::begin(thisTrackster.vertices()),
                    std::end(thisTrackster.vertices()),
                    std::back_inserter(outTrackster.vertices()));
          std::copy(std::begin(thisTrackster.vertex_multiplicity()),
                    std::end(thisTrackster.vertex_multiplicity()),
                    std::back_inserter(outTrackster.vertex_multiplicity()));
          candidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));
          if (chargedCommunities[i_c][i_cc]) {
            auto trackToTracksterEdge = graphs[i_c].getTrackToTracksterEdge();
            candidate.setTrackPtr(edm::Ptr<reco::Track>(tkH, trackToTracksterEdge.first));
          }
        }
        candidates.push_back(candidate);
        resultTrackstersMerged.push_back(outTrackster);
      }
      i_cc++;
    }
    i_c++;
  }

  assignPCAtoTracksters(
      resultTrackstersMerged, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());
  model.energyRegressionAndID(layerClusters, resultTrackstersMerged);
//  std::cout << " CLUE3D Trackster " << tracksters.size() << " Tracksters Merged " << resultTrackstersMerged.size()
            //<< std::endl;
}  // linkTracksters

void LinkingAlgoByLouvainAlgo::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("delta_ts_em_had", 0.03);
  desc.add<double>("delta_ts_had_had", 0.03);
  desc.add<double>("separationSmall", 2);  //cm
  desc.add<double>("separation", 6);       //cm
  desc.add<int>("maxDepth", 10);
  desc.add<double>("track_time_quality_threshold", 0.5);
  LinkingAlgoBase::fillPSetDescription(desc);
}
//loop over communities
/*    for (size_t i = 0; i != communities.size(); i++) {
      std::cout << "START COMMUNITY " << i << std::endl;
      Community &community = communities[i];
      float max_delta_modularity = 0;
      unsigned int best_community = i;
      std::cout << "Best community " << best_community << std::endl;
      //get node from communitiy
      for (size_t j = 0; j < community.nodes.size(); j++) {
        std::cout << "Cpmmunity node " << community.nodes[j] << "Trackster "
                  << graph.getNode(community.nodes[j]).getId() << std::endl;
        const Node &node_data = graph.getNode(community.nodes[j]);
        auto node_weight = node_data.getWeightedDegree();
        //loop over neighbour!
        for (auto &[neighbor, weight] : node_data.getWeightedEdges()) {
          int neighbor_community = node_to_community[trackster_to_node.at(neighbor)];
          if (neighbor_community == static_cast<int>(i)) {
            continue;
          }
          //now i have to remove the node from its community and put the node in the neighbouring community
          auto new_community_nodes = communities[neighbor_community].nodes;
          new_community_nodes.push_back(community.nodes[j]);

          std::cout << "Community internal weight " << community.internal_weight << " weight to add " << getWeightedEdgesInCommunity(node_data, i , node_to_community, trackster_to_node)
                    << std::endl;
          std::cout << "Community undergoing " << i << std::endl; 
          std::cout << "\t adding Node " << trackster_to_node.at(node_data.getId()) << " check internal weights " << std::endl;
          for(auto ed : node_data.getWeightedEdges()){
            std::cout << "\t\tcomm " << node_to_community.at(trackster_to_node.at(ed.first)) << std::endl;
          }
          float new_internal_weight = community.internal_weight + getWeightedEdgesInCommunity(node_data, i , node_to_community, trackster_to_node);
          community.internal_weight = new_internal_weight;
          std::cout << "Comm total weight " << community.total_weight << " add weight " << weight << " Weighted degree "
                    << graph.getNode(trackster_to_node.at(neighbor)).getWeightedDegree() << std::endl;
          float new_total_weight =
              community.total_weight +  2*graph.getNode(trackster_to_node.at(neighbor)).getWeightedDegree();
          std::cout << "new total weight " << new_total_weight << std::endl;
        community.total_weight = new_total_weight;
         // float new_modularity = new_internal_weight / new_total_weight -
           //                      std::pow(new_total_weight, 2) / std::pow(graph.getTotalWeight(), 2);
          std::cout << "Internal weight " << new_internal_weight << " total weight graph " << graph.getTotalWeight() << " comm total weight " << new_total_weight << std::endl;
          float new_modularity = (2 * new_internal_weight / graph.getTotalWeight()) -
                           std::pow(2 * new_total_weight / graph.getTotalWeight(), 2);

          std::cout << "New modularity " << new_modularity << " Modularity " << modularity(graph, node_to_community, trackster_to_node) << std::endl;
          std::cout << " Comm size " << communities.size() << " graph " << graph.size() << " NodeToComm " << node_to_community.size() << std::endl;
          float delta_modularity = new_modularity - modularity(graph, node_to_community, trackster_to_node);
          std::cout << "delta modularity " << delta_modularity << " max_delta_mod " << max_delta_modularity
                    << std::endl;
          if(community.internal_weight > community.total_weight){
            std::cout << "@@@@@@@@@@@@@@@@@@@ SOMETHING IS WRONG @@@@@@@@@@@@@@@" << std::endl;
            std::cout << "Internal " << community.internal_weight << " Total " << community.total_weight << std::endl;
          }
          if (delta_modularity > max_delta_modularity) {
            max_delta_modularity = delta_modularity;
            best_community = neighbor_community;
          }
          std::cout << " New best community " << best_community << std::endl;
        }
      }
      if (max_delta_modularity > 0) {
        changed = true;
        for (size_t j = 0; j < community.nodes.size(); j++) {
          const Node &node_data = graph.getNode(community.nodes[j]);
          std::cout << "Community Internal weight " << community.internal_weight << " to subratrct "
                    << getWeightedEdgesInCommunity(node_data, best_community, node_to_community, trackster_to_node)
                    << std::endl;

          community.internal_weight -=
              getWeightedEdgesInCommunity(node_data, best_community, node_to_community, trackster_to_node);
          std::cout << " Community total weight " << community.total_weight << " Node " << j << " degree "
                    << node_data.getWeightedDegree() << std::endl;
          community.total_weight -= node_data.getWeightedDegree();

          communities[best_community].internal_weight +=
              getWeightedEdgesInCommunity(node_data, best_community, node_to_community, trackster_to_node);
          communities[best_community].total_weight += node_data.getWeightedDegree();

          new_node_to_community[trackster_to_node.at(node_data.getId())] = best_community;
        }
      }
    }

    std::cout << "######################## END COMMUNITY " << std::endl;

    node_to_community = new_node_to_community;
    for (auto const &[node, comm] : node_to_community) {
      std::cout << "Node " << node << " Comm " << comm << std::endl;
    }
  
  */