#pragma once
#include <utility> //  std::swap
#include <cstddef> // std::size_t
#include <cmath> // std::sqrt
#include <vector>
#include "CommonTypes.h"
namespace BrainGraph {
  class Graph;
  /**
   * A property is mathematically a vector in some vector space defined on the
   * real numbers. 
   *
   * Should have:
   *   Implement arithmetic operators
   *
   * Nice to have:
   *   Add iterators so we can use range-based for loops. Iteration works fine
   *   Python, so it is of minor importance.
   *
   */
  class Property {
    friend class Graph;
  public:
    typedef property_type value_type;

    Property()
      : values()
    {}
    ~Property(){};

    Property(value_type value)
      : values(1, value)
    {}

    Property(std::vector<value_type> values)
      : values(values)
    {}

    Property(std::initializer_list<value_type> values)
      : values(values)
    {}

    Property(const Property& other)
    : values(other.values)
    {}

#ifndef SWIG
    // Stuff that swig does not like

    // swig does not like rvalues and consider it a syntax error, so we cant get a
    // move constructor. The assignment operator is ignored so there is no reason
    // to compile them or the swap function.
    // Swap resources
    friend void swap(Property& first, Property& second) {
      using std::swap;
      swap(first.values, second.values);
    }

    // Assignment constructor
    Property& operator=(Property other) {
      swap(*this, other);
      return *this;
    }

    // Move constructor
    Property(Property&& other)
    : Property() {
      swap(*this, other);
    }
#endif

    // Bounds checked element access
    value_type& operator[](std::size_t i) { return this->values.at(i); }
    // Const bounds checked element access
    const value_type& operator[](std::size_t i) const { return this->values.at(i); }
  
    // Dimension of the vector space
    std::size_t dim() const { return values.size(); }

    // Euclid norm
    value_type norm() const {
      value_type sum = 0;
      for (value_type value : values) {
	sum += value*value;
      }
      return std::sqrt(sum);
    }

    Property& operator-=(const Property& rhs) {
      for (std::size_t i = 0; i < this->dim(); ++i) {
	this->operator[](i) -= rhs[i];
      }
      return *this;
    }
    Property operator-(const Property& rhs) {
      Property lhs(*this);
      lhs -= rhs;
      return lhs;
    }

    Property& operator/=(value_type rhs) {
      for (std::size_t i = 0; i < this->dim(); ++i) {
	this->operator[](i) /= rhs;
      }
      return *this;
    }
    Property operator/(value_type rhs) {
      Property lhs(*this);
      lhs /= rhs;
      return lhs;
    }
    
    // In order to get seamless integration with python, it is necessary to define
    // the comparison operators as member functions. If we define the functions 
    // outside the class, we end up comparing pointer addresses or something
    // similar.
    // If someone figures out how to get comparison working when functions are
    // defined outside the class, feel free to change it.
    // 
    // An important thing about comparing vectors is that for two vectors v1, v2
    // !(v1 < v2) && !(v1 == v2) && !(v1 > v2)
    // could be true. This is because the norm is used for ordering, while
    // equality is component-wise
    bool operator==(const Property& rhs) const {
      if (this->dim() != rhs.dim()) {
	return false;
      }
      for (std::size_t i = 0; i < this->dim(); ++i) {
	if (this->operator[](i) != rhs[i]) {
	  return false;
	}
      }
      return true;
    }
    bool operator!=(const Property& rhs) const {return !(this->operator==(rhs));}

    // A Property is less than another property if the dimensions is lower,
    // or if the dimension is the same and the norm is lower.
    bool operator< (const Property& rhs) const { 
      return (this->dim() < rhs.dim() || 
	      (this->dim() == rhs.dim() && this->norm() < rhs.norm()));
    }
    bool operator> (const Property& rhs) const {return rhs.operator<(*this); }


    value_type dot(const Property& other) {
      value_type result = 0;
      for (std::size_t i = 0; i < this->dim(); ++i) {
	result += this->operator[](i) * other[i];
      }
      return result;
    }

  private:
    std::vector<value_type> values;
  };
}
