#pragma once
#include <iostream>
#include <functional>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <boost/heap/fibonacci_heap.hpp>
#include "CommonTypes.h"

namespace BrainGraph { 
  class Graph;
  class ROItoROI;

  /**
   * \brief Node type used in shortest path calculations.
   */
  struct ShortestPathNode {
    /**
     * ShortestPathNode constructor
     */ 
    ShortestPathNode(weight_type weight=0,
		     weight_type length=0,
		     id_type id=0,
		     id_type predecessor=0,
		     bool is_target=false) 
      : weight(weight),
	length(length),
	id(id),
	predecessor(predecessor),
	is_target(is_target) {}
    
    weight_type weight;  ///< Accumulated weight (=distance) from the source to this node.
    weight_type length;  ///< Number of edges from the source to this node.
    id_type id;          ///< Id of this node which matches the id in the graph.
    id_type predecessor; ///< Previous node on the path from source to this node.
    bool is_target;      ///< Flag indicating if this node is considered a required target.
  };
  typedef ShortestPathNode node_t;


  /**
   * \brief Comparison function used to implement a min-priority-queue on the weight of ShortestPathNodes.
   */
  struct CompareNode {
    bool operator()(const node_t* first, const node_t* second) const {
      return first->weight > second->weight;
    }
  };

  /*
   * Type definitions for the priority queue.
   */
  typedef boost::heap::compare<CompareNode> compare_t;
  typedef boost::heap::fibonacci_heap<node_t*, compare_t> priority_queue;
  typedef priority_queue::handle_type handle_t;

  /**
   * Calculates the single source shortest path, from source to AT LEAST all 
   * nodes that are considered targets. This may end up calculating all paths if
   * a target node is the farthest node.
   */ 
  std::vector<node_t> 
  single_source_shortest_path(const Graph &graph, 
			      id_type source, 
			      std::vector<bool> targets);

  /**
   * Returns the nodes that are in the path from the source used to calculate
   * the shortest paths, to the target node.
   * The path is returned in order from source to target,
   * such that if P is the path then
   * while (P[i].id != target) std::cout << P[i++].id;
   * will print the ids of the nodes going from source to target.
   */
  std::vector<node_t>
  path_to(std::vector<node_t> paths, 
	  id_type target);

  /**
   * Calculate shortest path between all nodes in sources and all nodes in targets.
   * The result is a ROItoROI structure. If r2r is the ROItoROI structure, then 
   * for all a_i in roi_a and all b_i in roi_b we have that
   * r2r[a_i][b_i]
   * is a vector of ShortestPathNodes containing the path from node a_i to node b_i.
   */
  ROItoROI 
  r2r_shortest_path(const Graph &graph,
		    std::vector<id_type> sources,
		    std::vector<id_type> targets);

}
