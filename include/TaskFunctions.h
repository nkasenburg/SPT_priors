#pragma once
#include <iostream>
#include <future>
#include <vector>
#include <chrono>

namespace TaskFunctions {
  /**
   * Wait untill all futures are ready
   * Implementation from:
   *   Bjarne Stroustrup : "The C++ Programming Language, Fourth Edition", page 1243.
   */ 
  template<typename T>
  std::vector<T>
  wait_for_all(std::vector< std::future<T> >& futures) {
    std::vector<T> results;
    for (auto& future : futures) {
      results.push_back(future.get());
    }
    return results;
  }

  /**
   * Wait untill any future in the range [begin, end) is ready
   * Implementation from
   *   Bjarne Stroustrup : "The C++ Programming Language, Fourth Edition", page 1243.
   */ 
  template<typename Iter>
  Iter 
  wait_for_any(Iter begin,
	       Iter end,
	       std::chrono::steady_clock::duration d) {
    while (true) {
      bool any_valid = false;
      for (auto it = begin; it != end; ++it) {
	if (it->valid()) {
	  any_valid = true;
	  switch(it->wait_for(std::chrono::seconds{0})) {
	  case std::future_status::ready:    return it;
	  case std::future_status::timeout:  break;
	  case std::future_status::deferred: throw std::runtime_error("wait_for_any(): deferred future");
	  }
	}
      }
      if (any_valid) {
	std::this_thread::sleep_for(d);
      } else {
	throw std::runtime_error("wait_for_any(): No valid task");
      }
    }
  }
}
