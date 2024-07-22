#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TICLGraph.h"
#include <unordered_set>
#include <vector>
#include <string>
#include <stack>

namespace ticl {

  void Node::findSubComponents(std::vector<Node>& graph, std::vector<unsigned int>& subComponent, std::string tabs) {
    std::stack<unsigned int> stack;
    stack.push(index_);

    while (!stack.empty()) {
      unsigned int currentIndex = stack.top();
      stack.pop();

      Node& currentNode = graph[currentIndex];
      if (!currentNode.alreadyVisited_) {
        LogDebug("TICLGraph") << tabs << " Visiting node " << currentNode.index_ << std::endl;
        currentNode.alreadyVisited_ = true;
        subComponent.push_back(currentIndex);

        for (auto const& neighbour : currentNode.outerNeighboursId_) {
          LogDebug("TICLGraph") << tabs << " Trying to visit " << neighbour << std::endl;
          if (!graph[neighbour].alreadyVisited_) {
            stack.push(neighbour);
          }
        }
      }
    }
  }
}//namespace ticl

std::vector<std::vector<unsigned int>> TICLGraph::findSubComponents() {
  std::vector<std::vector<unsigned int>> components;
  for (auto& node : nodes_) {
    if (!node.alreadyVisited()) {
      if (node.hasOuterNeighbours() || (!node.hasOuterNeighbours() && !node.hasInnerNeighbours())) {
        std::string tabs = "\t";
        std::vector<unsigned int> tmpSubComponents;
        node.findSubComponents(nodes_, tmpSubComponents, tabs);
        components.push_back(tmpSubComponents);
      }
    }
  }
  return components;
}
void TICLGraph::dfsForCC(unsigned int nodeIndex,
                         std::unordered_set<unsigned int>& visited,
                         std::vector<unsigned int>& component) const {
  visited.insert(nodeIndex);
  component.push_back(nodeIndex);

  for (auto const& neighbourIndex : nodes_[nodeIndex].getOuterNeighbours()) {
    if (visited.find(neighbourIndex) == visited.end()) {
      dfsForCC(neighbourIndex, visited, component);
    }
  }
}

std::vector<std::vector<unsigned int>> TICLGraph::getConnectedComponents() const {
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
