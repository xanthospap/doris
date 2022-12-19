#include "var_utils.hpp"

int dso::get_yaml_value_depth2(const YAML::Node &root, const char *key1,
                               const char *key2, char *buf) noexcept {
  try {
    const YAML::Node node = dso::get_yaml_node(root, key1);
    std::string result = node[std::string(key2)].as<std::string>();
    std::strcpy(buf, result.c_str());
  } catch (...) {
    fprintf(stderr, "ERROR@%s Failed finding key %s::%s\n", __func__, key1,
            key2);
    return 1;
  }
  return 0;
}
