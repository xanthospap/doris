#include "icgemio.hpp"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

constexpr int max_header_lines = 1000;
constexpr int bsz = 256;

#ifdef DEBUG
void dso::Icgem::print_details() {
  printf("Details of igem parsed header for file: %s\n", filename.c_str());
  printf("Product Type : %s\n", product_type);
  printf("Model Name   : %s\n", modelname);
  printf("Tide System  : %s\n", tide_system);
  printf("Normalization: %s\n", norm);
  printf("Errors       : %s\n", errors);
  printf("GM           : %15.10f\n", earth_gravity_constant);
  printf("R            : %15.10f\n", radius);
  printf("Max Degree   : %d\n", max_degree);
  return;
}
#endif

/// Given a C-String, go to the first non-whitespace character. If the line is
/// empty or only contains whitespace characters, the returned (pointer to
/// char) is NULL
/// example: given line="  foo bar\0"
///          return a pointer to 'f', aka the string "foo bar"
/// example: given "   \0"
///          return '\0'
const char *next_non_ws_char(const char *line) noexcept {
  const char *s = line;
  while (*s && *s == ' ')
    ++s;
  return s;
}

const char *next_ws_char(const char *line) noexcept {
  const char *s = line;
  while (*s && *s != ' ')
    ++s;
  return s;
}

/// Check if a given line corresponds to a given parameter name (or keyword)
/// parameter_name Keyword to match (e.g. 'product_type', 'modelname', etc)
/// line           Line to check
/// keyword_size
const char *resolve_header_keyword(const char *parameter_name, const char *line,
                                   int &keyword_size) noexcept {
  auto psz = std::strlen(parameter_name);
  keyword_size = 0;

  // does the line start with the parameter_name string ?
  if (!std::strncmp(line, parameter_name, psz)) {
    // start of new word after the keyword in line
    const char *start = next_non_ws_char(line + psz + 1);
    // start of second word after the keyword in line
    const char *stop = next_ws_char(start);
    if (!start || !stop) {
      keyword_size = -1;
      return nullptr;
    }
    keyword_size = stop - start + 1;
    return start;
  } else {
    // parameter_name not matched in begining of line; return NULL
    return nullptr;
  }
}

/// @warning this function assumes that comment and header lines do not exceed
///          256 characters; this is not guaranteed by the icgem format, but
///          practically holds for all cases.
int dso::Icgem::parse_header() noexcept {
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening icgem file %s (traceback: %s)\n",
            filename.c_str(), __func__);
    return 1;
  }

  char buf[bsz];
  const char *keyword;
  int keyword_size;
  int counter = 0;
  int error = 0;

  // keep reading new lines untill we meet a line starting with 'end_of_head'
  fin.getline(buf, bsz);
  while (std::strncmp(buf, "end_of_head", 11)) {
    // the line read can start with leading whitespace characters; skip them
    // it can also be empty!
    const char *line = next_non_ws_char(buf);

    // check if this the 'product_type' field (mandatory keyword)
    keyword = resolve_header_keyword("product_type", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "product_type", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      assert(keyword_size < (int)sizeof(product_type));
      std::strncpy(product_type, keyword, keyword_size);
    }

    // check if this the 'modelname' field (mandatory keyword)
    keyword = resolve_header_keyword("modelname", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "modelname", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      assert(keyword_size < (int)sizeof(modelname));
      std::strncpy(modelname, keyword, keyword_size);
    }

    // check if this the 'earth_gravity_constant' field (mandatory keyword)
    keyword =
        resolve_header_keyword("earth_gravity_constant", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "earth_gravity_constant", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      char *end;
      earth_gravity_constant = std::strtod(keyword, &end);
      if (end == keyword || earth_gravity_constant == 0e0) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "earth_gravity_constant", filename.c_str(), __func__);
        error = 2;
      }
    }

    // check if this the 'radius' field (mandatory keyword)
    keyword = resolve_header_keyword("radius", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "radius", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      char *end;
      radius = std::strtod(keyword, &end);
      if (end == keyword || radius == 0e0) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "radius", filename.c_str(), __func__);
        error = 2;
      }
    }

    // check if this the 'max_degree' field (mandatory keyword)
    keyword = resolve_header_keyword("max_degree", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "max_degree", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      char *end;
      max_degree = std::strtol(keyword, &end, 10);
      if (end == keyword || max_degree == 0e0) {
        fprintf(stderr,
                "[ERROR] Failed parsing value for parameter %s in icgem "
                "%s(traceback: %s)\n",
                "max_degree", filename.c_str(), __func__);
        error = 2;
      }
    }

    // check if this the 'errors' field (mandatory keyword)
    keyword = resolve_header_keyword("errors", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "errors", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      assert(keyword_size < (int)sizeof(errors));
      std::strncpy(errors, keyword, keyword_size);
    }

    // check if this the 'tide_system' field (optional keyword)
    keyword = resolve_header_keyword("tide_system", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "tide_system", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      assert(keyword_size < (int)sizeof(tide_system));
      std::strncpy(tide_system, keyword, keyword_size);
    }

    // check if this the 'norm' field (optional keyword)
    keyword = resolve_header_keyword("norm", line, keyword_size);
    if (keyword_size < 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing value for parameter %s in icgem "
              "%s(traceback: %s)\n",
              "norm", filename.c_str(), __func__);
      error = 1;
    } else if (keyword) {
      assert(keyword_size < (int)sizeof(norm));
      std::strncpy(norm, keyword, keyword_size);
    }

    fin.getline(buf, bsz);
    if (error || (++counter > max_header_lines)) {
      if (error)
        return error;
      fprintf(stderr,
              "[ERROR] Failed parsing header in icgem %s; too many header "
              "lines ... (traceback: %s)\n",
              filename.c_str(), __func__);
      return 1;
    }

  } // keep on parsing (header) lines ...

  data_section_pos = fin.tellg();
  return 0;
}
