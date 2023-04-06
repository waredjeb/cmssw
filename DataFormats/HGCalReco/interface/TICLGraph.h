#ifndef RecoHGCal_TICLGraph_h
#define RecoHGCal_TICLGraph_h

#include <vector>
#include <utility>

class Node {
public:
  Node() = default;
  Node(unsigned index, bool isTrackster = true) : index_(index), isTrackster_(isTrackster){};

  void addInner(unsigned int trackster_id) { innerNodes_.push_back(trackster_id); }
  void addOuter(unsigned int trackster_id) { outerNodes_.push_back(trackster_id); }
  void addEdge(unsigned int trackster_id) { edges_.push_back(trackster_id); }
  void addEdge(unsigned int trackster_id, const float weight) { edgesWeight_.emplace_back(trackster_id, weight); }
  const unsigned int getId() const { return index_; }
  bool isTrackster() const { return isTrackster_; }
  std::vector<unsigned int> getInner() const { return innerNodes_; }
  std::vector<unsigned int> getOuter() const { return outerNodes_; }
  std::vector<unsigned int> getEdges() const { return edges_; }
  std::vector<std::pair<unsigned int, float>> getWeightedEdges() const { return edgesWeight_; }
  float getWeightedDegree() const {
      float degree = 0;
      for (const auto& [neighbor, weight] : edgesWeight_) {
          //degree += weight;
          degree += weight;

      }
      return degree;
  }

  ~Node() = default;

private:
  unsigned index_;
  bool isTrackster_;
  std::vector<unsigned int> innerNodes_;
  std::vector<unsigned int> outerNodes_;
  std::vector<unsigned int> edges_;
  std::vector<std::pair<unsigned int, float>> edgesWeight_;
};

class TICLGraph {
public:
  TICLGraph() = default;
  TICLGraph(std::vector<Node> &n) :
    nodes_ {n}
  {

    for(auto const& n : nodes_){
      for (auto const& [e,w] : n.getWeightedEdges()){
        totalWeight_+= w;
      }
    }
  };

  void setTrackToTracksterEdge(const int track_index, const int trackster_index){
    track_to_trackster_.first = track_index;
    track_to_trackster_.second = trackster_index;
  }
  const std::pair<int,int> getTrackToTracksterEdge() const {
    return track_to_trackster_;
  }
  const std::vector<Node>& getNodes() const { return nodes_; }
  const Node& getNode(unsigned int i) const { return nodes_[i]; }
  void addNode(const Node& node) { nodes_.push_back(node); }
  unsigned int size() const { return nodes_.size(); }
  float getTotalWeight() const { return totalWeight_; }

  std::vector<std::pair<unsigned int, float>> getNeighborsWithWeights(unsigned int i) const {
    return nodes_[i].getWeightedEdges();
  }

  float getDegree(unsigned int i) const {
    return nodes_[i].getWeightedDegree();
  }

  float getSelfWeight(unsigned int i) const {
    for (const auto& [neighbor, weight] : nodes_[i].getWeightedEdges()) {
      if (neighbor == i) {
        return weight;
      }
    }
    return 0.0;
  }

  float getEdgeWeight(unsigned int i, unsigned int j) const {
  for (const auto& [neighbor, weight] : nodes_[i].getWeightedEdges()) {
    if (neighbor == j) {
      return weight;
    }
  }
  return 0.0;
}


  ~TICLGraph() = default;

private:
  std::vector<Node> nodes_;
  std::pair<int, int> track_to_trackster_ = std::make_pair<int,int>(-1,-1); 
  float totalWeight_ = 0.0;

};

#endif
