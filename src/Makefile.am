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
-std=c++17 \
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
  sp3.hpp

##
##  Source files (distributed).
## ----------------------------
##
dist_libcppdoris_la_SOURCES = \
  doris_system_info.cpp \
  read_rnx_header.cpp \
  doris_rinex.cpp \
  sp3_read_header.cpp \
  sp3.cpp
