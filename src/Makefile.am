##
##  The library will be called: (lib)"cppdoris".
## ----------------------------------------------
##
lib_LTLIBRARIES = libcppdoris.la

##
##  Only use the following flags for debuging purposes.
## ----------------------------------------------------
##
libcppdoris_la_CXXFLAGS = \
  -std=c++20 \
  -g \
  -pg \
  -Wall \
  -Wextra \
  -Werror \
  -pedantic \
  -W \
  -Wshadow \
  -Winline \
  -Wdisabled-optimization \
  -march=native \
  -DDEBUG

##
##  Header files (distributed) installed in /$(includedir)/$(package name).
## ------------------------------------------------------------------------
##
pkginclude_HEADERS = \
  doris_system_info.hpp \
  doris_rinex.hpp \
  antenna_pcv.hpp \
  doris_utils.hpp

##
##  Source files (distributed).
## ----------------------------
##
dist_libcppdoris_la_SOURCES = \
  doris_system_info.cpp \
  read_rnx_header.cpp \
  doris_rinex.cpp \
  extract_doppler_count.cpp \
  get_pdop_crd.cpp
