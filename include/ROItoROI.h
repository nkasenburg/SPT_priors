#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "Dijkstra.h"

namespace BrainGraph {
  /**
   * A wrapper aroud a 2-level deep map (aka map<key, map<key, value>> ) that
   * supports iterating over the values in the inner map in lexicographic 
   * order. So if we use int as key, then
   *
   * for (int i : sort(outer_keys))
   *   for (int j : sort(inner_keys[i]))
   *     // do something with r2r(i,j)
   *
   * is equivalent to
   *
   * for (auto value : r2r)
   *  // do something with the value
   */

//template<typename ROIType, typename PathType>
  typedef unsigned long ROIType;
  typedef std::vector<BrainGraph::ShortestPathNode> PathType;
  class ROItoROI {
  public:
    // SWIG does not support the using method for aliasing
    typedef std::map<ROIType, PathType> InnerMap;
    typedef std::map<ROIType, InnerMap> OuterMap;

    typedef ROIType key_type;
    typedef PathType value_type;
    typedef PathType& reference;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t size_type;

    // using InnerMap = std::map<unsigned long, std::vector<ShortestPathNode> >;
    // using OuterMap = std::map<unsigned long, InnerMap>;
  
    // using value_type = std::vector<ShortestPathNode> ;
    // using reference = std::vector<ShortestPathNode> &;
    // using iterator_category = std::forward_iterator_tag;
    // using size_type = typename OuterMap::size_type;

    
    class iterator {
      friend ROItoROI;
    public:
      iterator() {}
      iterator(const iterator& other) 
        : outer(other.outer)
	, inner(other.inner)
	, inner_end(other.inner_end)
      {}

      ~iterator() {}

      bool operator==(const iterator& other) const { 
	return this->outer == other.outer && 
	  (this->inner == other.inner || this->inner == this->inner_end);
      }
      bool operator!=(const iterator& other) const { return !this->operator==(other); }

      iterator& operator++() {
	++(this->inner);
	if (this->inner == this->inner_end) {
	  ++(this->outer);
	  if (this->outer != this->outer_end) {
	    this->inner = this->outer->second.begin();
	    this->inner_end = this->outer->second.end();
	  }
	}
	return *this;
      }

      reference operator*() { return this->inner->second; }
      
    private:
      typename OuterMap::iterator outer;
      typename InnerMap::iterator inner;
      typename OuterMap::iterator outer_end;
      typename InnerMap::iterator inner_end;

      iterator(ROItoROI& r2r, bool end) {
	if (end) {
	  this->outer = this->outer_end = r2r.paths.end();
	  this->inner = this->inner_end = r2r.paths.begin()->second.end();
	} else {
	  this->outer = r2r.paths.begin();
	  this->outer_end = r2r.paths.end();
	  this->inner = outer->second.begin();
	  this->inner_end = outer->second.end();
	}
      }
    };

    ROItoROI() : paths() {};
    // ROItoROI(std::initalizer_list) : paths() {};
    ~ROItoROI() {};

    // const std::vector<ShortestPathNode>& operator()(key_type from, key_type to) const { 
    //   return this->paths[from][to];
    // }
    reference operator()(key_type from, key_type to) { 
      try {
	return this->paths.at(from)[to];
      } catch (std::out_of_range) {
	this->paths[from] = InnerMap();
	return this->paths[from][to];
      }
    }

    InnerMap& operator()(key_type from) { return this->paths[from]; }
    //const InnerMap& operator()(key_type from) const { return this->paths[from]; }

    reference get_path(key_type from, key_type to) { 
      return this->operator()(from,to);
    }
    void set_path(key_type from, key_type to, value_type path) { 
      this->operator()(from,to) = path;
    }

    //InnerMap& get_inner(key_type from) { return this->operator()(from); }

    size_type no_of_sources() const { return this->paths.size(); }
    size_type no_of_targets(key_type k) const { return this->paths.at(k).size(); }
    
    iterator begin() { return iterator(*this, false); }
    iterator end() { return iterator(*this, true); }

  private:
    OuterMap paths;
  };
}
