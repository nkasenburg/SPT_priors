#pragma once
#include "CommonTypes.h"

namespace BrainGraph {
  /**
   * An edge containing one end of the edge and the weight of the edge. 
   */ 
  struct Edge {
    Edge()
      : node(0)
      , weight(0)
    {}

    Edge(id_type node, weight_type weight)
      : node(node)
      , weight(weight)
    {}

    bool operator==(const Edge& rhs) const {
      return this->node == rhs.node && this->weight == rhs.weight;
    }

    bool operator!=(const Edge& rhs) const { return !this->operator==(rhs);}

    id_type node;       ///< One end of the edge
    weight_type weight; ///< The weight of the edge
  };
}
