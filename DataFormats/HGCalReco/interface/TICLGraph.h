#ifndef DataFormats_HGCalReco_TICLGraph_h
#define DataFormats_HGCalReco_TICLGraph_h

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <unordered_set>

class Node {
public:
  Node() = default;
  Node(unsigned index, bool isTrackster = true) : index_(index), isTrackster_(isTrackster){};

  void addNeighbour(unsigned int trackster_id) { neighbours_.push_back(trackster_id); }

  const unsigned int getId() const { return index_; }
  std::vector<unsigned int> getNeighbours() const { return neighbours_; }

    
  ~Node() = default;

private:
  unsigned index_;
  bool isTrackster_;
  std::vector<unsigned int> neighbours_;
};

class TICLGraph {
public:
 // can i remove default constructor ?? edm::Wrapper problem
 // without default constructor i could initialize connectedComponents when building the Graph
  TICLGraph() = default;
  TICLGraph(std::vector<Node> &n) {
      nodes_ = n;
  };
  const std::vector<Node> &getNodes() const { return nodes_; }
  const Node &getNode(int i) const { return nodes_[i]; }
  ~TICLGraph() = default;

  void dfsForCC(unsigned int nodeIndex, std::unordered_set<unsigned int>& visited, std::vector<unsigned int>& component) const{
    visited.insert(nodeIndex);
    component.push_back(nodeIndex);

    for(auto const& neighbourIndex : nodes_[nodeIndex].getNeighbours()){
        if(visited.find(neighbourIndex) == visited.end()){
            dfsForCC(neighbourIndex, visited, component);
        }
    }
  }

  std::vector<std::vector<unsigned int>> getConnectedComponents() const{
      std::unordered_set<unsigned int> visited;
      std::vector<std::vector<unsigned int>> components;
    
        for (unsigned int i = 0; i < nodes_.size(); ++i) {
            if (visited.find(i) == visited.end()) {
                std::vector<unsigned int> component;
                dfsForCC(i, visited, component);
                components.push_back(component);
            }
        }
        
        return components;
    }
    
private:
  std::vector<Node> nodes_;
};

#endif
