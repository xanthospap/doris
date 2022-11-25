#ifndef __DSO_DORIS_VAR_UTILS_HPP__
#define __DSO_DORIS_VAR_UTILS_HPP__

#include "var/celestrak.hpp"
#include "yaml-cpp/yaml.h"
#include <cstring>
#include <exception>

namespace dso {

/// @brief Given the parent node root, fetch a child node
/// @param[in] root The parent node (YAML::Node)
/// @param[in] key  The key of the child node to return
/// @return a YAML::Node, child of the root node, specified by key "key"
inline YAML::Node get_yaml_node(const YAML::Node &root, const char *key) {
  return root.operator[](std::string(key));
}

/// @brief Get the value of a depth-2 'dictionary' key/value pair
/// @param[in] root The parent node (YAML::Node)
/// @param[in] key1  The key of the child node to return
/// @param[in] key2  The key for the value we want
/// @param[out] val The resolved value
/// @return Anything other than 0 denotes an error and val should not be used.
/// Example (in some yaml file 'conf.yaml')
/// ---
/// foo:
///  bar1: 1
///  bar2: 2
/// data:
///   a-coeff: 3
///   b-coeff: 4
/// ...
/// const YAML::Node config = YAML::LoadFile(conf.yaml);
/// if (get_yaml_value_depth2(config, "data", "b-coeff", val))
///   // handle error ...
/// else
///   // val hold the value 4
template <typename T>
int get_yaml_value_depth2(const YAML::Node &root, const char *key1,
                          const char *key2, T &val) noexcept {
  try {
    const YAML::Node node = get_yaml_node(root, key1);
    val = node[std::string(key2)].as<T>();
  } catch (...) {
    fprintf(stderr, "ERROR@%s Failed finding key %s::%s\n", __func__, key1,
            key2);
    return 1;
  }
  return 0;
}

/// @brief Get the value of a depth-2 'dictionary' key/value pair as C-string
/// @param[in] root The parent node (YAML::Node)
/// @param[in] key1  The key of the child node to return
/// @param[in] key2  The key for the value we want
/// @param[out] buf  The resolved value as a C-string
/// @return Anything other than 0 denotes an error
/// Example (in some yaml file 'conf.yaml')
/// ---
/// foo:
///  bar1: 1
///  bar2: 2
/// data:
///   file1: foo/bar
///   file2: bar/foo
/// ...
/// const YAML::Node config = YAML::LoadFile(conf.yaml);
/// if (get_yaml_value_depth2(config, "data", "b-coeff", buf))
///   // handle error ...
/// else
///   // buf hold the string "bar/foo"
int get_yaml_value_depth2(const YAML::Node &root, const char *key1,
                          const char *key2, char *buf) noexcept;

/// @brief Get the value of a depth-3 'dictionary' key/value pair
/// @param[in] root The parent node (YAML::Node)
/// @param[in] key1  The key1 of the first child node (depth 1)
/// @param[in] key1  The key2 of the second child node (depth 2)
/// @param[in] key2  The key3 for the value we want
/// @param[out] val The resolved value
/// @return Anything other than 0 denotes an error and val should not be used.
/// Example (in some yaml file 'conf.yaml')
/// ---
/// foo:
///  bar1: 1
///  bar2: 2
/// data:
///   a-coeff: 3
///   b-coeff: 4
/// model:
///   foo:
///   bar:
///     koko: 1
/// ...
/// const YAML::Node config = YAML::LoadFile(conf.yaml);
/// if (get_yaml_value_depth3(config, "model", "bar", "koko", val))
///   // handle error ...
/// else
///   // val hold the value 1
template <typename T>
int get_yaml_value_depth3(const YAML::Node &root, const char *key1,
                          const char *key2, const char *key3, T &val) {
  const YAML::Node node = get_yaml_node(root, key1);
  return get_yaml_value_depth2(node, key2, key3, val);
}

/// @brief Specialization of the template <typename T>::get_yaml_value_depth3
/// function for C-strings
inline int get_yaml_value_depth3(const YAML::Node &root, const char *key1,
                                 const char *key2, const char *key3,
                                 char *buf) {
  const YAML::Node node = get_yaml_node(root, key1);
  return get_yaml_value_depth2(node, key2, key3, buf);
}

// int get_CelesTrack_flux_data(
//     const dso::modified_julian_day mjd, const char *fncsv,
//     dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {
//   return dso::utils::celestrak::details::parse_csv_for_date(mjd, fncsv,
//                                                              flux_data);
// }
//
// template <typename T>
// int get_CelesTrack_flux_data(
//     const dso::datetime<T> &t, const char *fncsv,
//     dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {
//   return dso::utils::celestrak::details::parse_csv_for_date(t.as_mjd(),
//   fncsv,
//                                                              flux_data);
// }

} // namespace dso

#endif
