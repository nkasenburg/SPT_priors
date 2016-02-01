#pragma once
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include "Property.h"
#include "CommonTypes.h"

namespace BrainGraph {
  class Graph;

  /**
   * A node with an id and a map of properties
   */
  class Node {
    friend class Graph;
  public:
    Node()
      : id_()
      , properties()
    {}

    ~Node(){}

    Node(id_type id)
      : id_(id)
      , properties()
    {}

    Node(id_type id, std::map<std::string, Property> properties)
      : id_(id)
      , properties(properties)
    {}
    
    Node(id_type id, std::initializer_list< std::map<std::string, Property>::value_type > properties)
      : id_(id)
      , properties(properties)
    {}

    Node(const Node& other) 
    : id_(other.id())
    , properties(other.properties)
    {}
    
    
#ifndef SWIG
    friend void swap(Node& first, Node& second) {
      using std::swap;
      swap(first.id_, second.id_);
      swap(first.properties, second.properties);
    }
    
    Node(Node&& other)
    : Node(0) {
      swap(*this, other);
    }
  
    Node& operator=(Node other) {
      swap(*this, other);
      return *this;
    }
#endif

    id_type id() const { return this->id_; }

    // In python equality et al are member functions, so to make the python
    // interface nicer, they are defined as member functions.
    bool operator==(const Node& rhs) const {
      if (this->id() != rhs.id() || this->properties.size() != rhs.properties.size()) {
	return false;
      }
      for (const auto& kv : this->properties) {
	try {
	  if (kv.second != rhs.properties.at(kv.first)) {
	    return false;
	  }
	} catch (std::out_of_range) {
	  return false;
	}
      }
      return true;
    }
    bool operator!=(const Node& rhs) const { return !this->operator==(rhs);}

  private:
    id_type id_;
  public:
    std::map<std::string, Property> properties;
  };
}
