#include <stdexcept>
#include <iostream>
#include <limits> // infinity()
#include <cmath> // M_PI
#include <functional>
#include <cassert>
#include <algorithm>
#include "Dijkstra.h"
#include "GeneralFunctions.h"
#include "Graph.h"

namespace BrainGraph { 

  // Single source shortest path (Dijkstra), stops when all target nodes have been reached.
  std::vector<node_t>
  single_source_shortest_path(const Graph &graph, 
			      id_type source, 
			      std::vector<bool> targets) {
    std::size_t target_count = 0;
    priority_queue unvisited;
    // We must reserve space for the nodes, otherwise the nodes will be moved when
    // the vector is resized and the pointers in the the priority_queue will be
    // invalid.
    std::vector<node_t> nodes(graph.no_of_nodes()); 
    std::vector<handle_t> handles(graph.no_of_nodes());
    
    for (id_type i = 0; i < graph.no_of_nodes(); ++i) {
      if (targets[i]) { ++target_count; }
      nodes[i] = {std::numeric_limits<weight_type>::infinity(), 0, i, i, targets[i]};
      handles[i] = unvisited.push(&nodes[i]);
    }
    nodes[source].weight = 0;
    unvisited.increase(handles[source]);
    while (! unvisited.empty()) {
      auto* current = unvisited.top();
      if (current->is_target && --target_count == 0) { break; } // Break if this is the last one
      unvisited.pop();
      for (auto edge : graph.edges(current->id)) {
  	auto& neighbour = nodes[edge.node];
  	if (current->weight + edge.weight < neighbour.weight) {
  	  neighbour.weight = current->weight + edge.weight;
  	  unvisited.increase(handles[neighbour.id]); // We decrease the weight, but increase the priority
  	  neighbour.length = current->length + 1;
	  neighbour.predecessor = current->id;
  	}
      }
    }
    return nodes;
  }


  // Get the path from source to target based on Dijkstra output
  std::vector<node_t> path_to(std::vector<node_t> paths, id_type target) {
    std::vector<node_t> this_path(paths[target].length+1);
    auto current = target;
    std::size_t i = this_path.size() - 1;
    //iterate through predecessors (starting from the target) until reaching the source
    while (paths[current].predecessor != current) {
      this_path[i--] = paths[current];
      current = paths[current].predecessor;
    }
    assert(i == 0);
    this_path[i] = paths[current];
    return this_path; 
  }

  // Get all shortest paths from a single source to many targets
  std::pair<id_type, ROItoROI::InnerMap>
  single_source_many_targets(const Graph &graph,
			     id_type source, 
			     std::vector<id_type> targets) {
    ROItoROI::InnerMap paths;
    //create target mask
	std::vector<bool> target_mask(graph.no_of_nodes());
	assert(std::none_of(target_mask.begin(), target_mask.end(),
    			GeneralFunctions::id<bool>)); // Verify all are false
    for (auto target : targets) { target_mask[target] = true; }
    //run Dijsktra for source
    auto shortest_paths = single_source_shortest_path(graph, source, target_mask);
    //iterate through targets and construct the corresponding paths
    for (auto target : targets) {
    	paths[target] = path_to(shortest_paths, target);
    }
    return std::make_pair(source, paths);
  }

  // Get all shortest paths from sources to targets
  ROItoROI 
  r2r_shortest_path(const Graph &graph, 
		    std::vector<id_type> sources, 
		    std::vector<id_type> targets) {
    ROItoROI paths;
    //iterate through all source voxels and compute the shortest paths to all target voxels
    for (auto source : sources) {
    	paths(source) = single_source_many_targets(graph, source, targets).second;
    }
    return paths;
  }

}
