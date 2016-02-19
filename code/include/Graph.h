#pragma once
#include <map>
#include <stdexcept>
#include <vector>
#include <utility>
#include "Node.h"
#include "Property.h"
#include "Edge.h"
#include "ROItoROI.h"
#include "CommonTypes.h"

namespace BrainGraph {
  /**
   * A graph class.
   * Supports both directed and undirected edges.
   * 
   */
  class Graph {
  public:

    /**
     * Load a graph stored in binary format.
     *
     * \param path  Path to graph file
     * \return      Graph
     */ 
    static Graph
    load_binary(std::string path);

    /**
     * Construct an empty graph
     */ 
    Graph()
      : nodes_()
      , edges_()
    {}

    /**
     * Copy construct a graph.
     * Nodes and edges are copied, so if they do not contain pointers/references
     * there will be no shared state.
     */ 
    Graph(const Graph& other)
    : nodes_(other.nodes_)      
    , edges_(other.edges_)
    {}

    /**
     * Graph constructor
     * Construct a graph from a node vector and an edge vector. The size of the
     * two vectors should be equal.
     *
     * For all nodes, the id of the node is changed so all nodes are assigned a
     * unique id in the range [0;nodes.size()]. This new id is referred to as
     * the index or idx. The edges are updated so they refer to this new index.
     *
     * It is possible to retrieve a node from the graph using either the
     * original id or the new index. Using the original id is in general slower
     * than using the new index.
     */
    Graph(std::vector<Node> nodes, std::vector< std::vector<Edge> > edges)
      : nodes_(nodes)
      ,	edges_(edges)
    {
      if (edges.size() != nodes.size()) {
	throw std::invalid_argument("There must be exactly one edgelist for each node. Got " + 
				    std::to_string(edges.size()) + " edges, and " +
				    std::to_string(nodes.size()) + " nodes.");
      }
    }

    
    /**
     * Graph destructor
     */
    ~Graph(){}


    /**
     * Get a node by the index assigned to the node
     * Time complexity: O(1)
     *
     * \param idx  Id of the node to get
     */
    const Node&
    node(id_type idx) const { return this->nodes_.at(idx); }
    
    Node& 
    node(id_type idx) { return this->nodes_.at(idx); }


    /**
     * Get the outgoing edges from a node identified by the index assigned to the node
     * Time complexity: O(log(number of nodes))
     * 
     * \param idx  Id of node to get edges from
     */
    const std::vector<Edge>&
    edges(id_type idx) const { return this->edges_.at(idx); }
    
    std::vector<Edge>& 
    edges(id_type idx) { return this->edges_.at(idx); }


    /**
     * Get the edge going from one node to another. Both nodes are identified by
     * the index assigned to them.
     * Time complexity: O(number of edges from node)
     */
    Edge&
    edge(id_type from_idx, id_type to_idx) {
      for (auto& edge: this->edges_.at(from_idx)) {
    	if (edge.node == to_idx) {
    	  return edge;
    	}
      }
      throw std::out_of_range("No edge from " + std::to_string(from_idx) + " to " + std::to_string(to_idx));
    }

    const Edge&
    edge(id_type from_idx, id_type to_idx) const {
      for (auto& edge: this->edges_.at(from_idx)) {
    	if (edge.node == to_idx) {
    	  return edge;
    	}
      }
      throw std::out_of_range("No edge from " + std::to_string(from_idx) + " to " + std::to_string(to_idx));
    }


    /**
     * Get the number of nodes in the graph
     * Time complexity: O(1)
     */
    id_type
    no_of_nodes() const { return nodes_.size(); }


    /**
     * Save the graph to the given path in binary format.
     *
     * \param path  Location where graph should be stored
     */
    bool
    save_binary(std::string path) const;

    /**
     * Calculate the importance of the node in the ROItoROI.
     * The number of times a node is used in a path in the ROItoROI is
     * node.properties["count"]
     * The relative weight of a node is stored in
     * node.properties["confidence"]
     */ 
    void
    calculate_node_importance(BrainGraph::ROItoROI r2r,    ///< ROItoROI used to calculate importance
			      std::string property_to_use_as_weight="prior",           ///< Key of node property to use as prior
			      std::string property_to_store_count_in="count",          ///< Key of node property to store count in
			      std::string property_to_store_importance_in="confidence" ///< Key of node property to store confidence in
			      );

    /**
     * Uses nodes position property to arrange nodes in a 3D matrix.
     * The values at each point in the matrix is given by node.properties[property_to_use_as_weight].
     * The matrix is returned in row-major order
     * x is the row dimension, y is column and z is slice
     * index(x,y,z) = z + y * slice + x * column * slice
     * position(i)  = (x,y,z) 
     *    where x   = i / (column * slice)
     *          y   = rem / slice
     *          z   = rem % slice
     *          rem = i % (column * slice)
     */
    std::vector<property_type>
    as_matrix(id_type rows,    ///< Number of rows in matrix
	      id_type columns, ///< Number of columns in matrix
	      id_type slices,  ///< Number of slices in matrix
	      std::string property_to_use_as_weight="confidence",      ///< Key of node property to use as confidence
	      std::string property_to_use_as_position="position" ///< Key of node property to use as position
	      ) const;
    
  private:
    std::vector<Node> nodes_;
    std::vector< std::vector<Edge> > edges_;
  };
}

