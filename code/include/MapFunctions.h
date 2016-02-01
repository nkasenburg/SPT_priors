#pragma once
#include <map>
#include <vector>
#include <algorithm>
#include <utility>

namespace MapFunctions {
  /**
   * Get the keys that are used in the map
   */
  template<typename Key, typename Value>
  std::vector<Key>
  keys(std::map<Key, Value> map) {
    std::vector<Key> keys(map.size());
    auto get_key = [](std::pair<const Key, Value> kv) -> Key { return kv.first; };
    std::transform(map.begin(), map.end(), keys.begin(), get_key);
    return keys;
  }

  /**
   * Get the values that are used in the map
   */
  template<typename Key, typename Value>
  std::vector<Value>
  values(std::map<Key, Value> map) {
    std::vector<Value> values(map.size());
    auto get_value = [](std::pair<const Key, Value> kv) -> Value { return kv.second; };
    std::transform(map.begin(), map.end(), values.begin(), get_value);
    return values;
  }

  /**
   * Check if a key is in a map
   */
  template<typename Key, typename Value>
  bool
  member(Key key, std::map<Key, Value> map) {
    return map.find(key) != map.end();
  }

  
  /**
   * Get keys in map sorted by their mapped value
   */
  template<typename Key, typename Value>
  std::vector<Key>
  sort_by_value(std::map<Key, Value> map) {
    using value_type = typename std::map<Key, Value>::value_type;

    std::vector< std::pair<Value, Key> > vk(map.size());
    auto T = [](value_type kv) { return std::make_pair(kv.second, kv.first); };
    std::transform(map.begin(), map.end(), vk.begin(), T);
    std::sort(vk.begin(), vk.end());

    std::vector<Key> keys(vk.size());
    std::transform(vk.begin(), vk.end(), keys.begin(), [](std::pair<Value, Key> x) { return x.second; });

    return keys;
  }
}
