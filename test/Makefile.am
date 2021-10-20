noinst_PROGRAMS = testDorisRinex testAntennaPcv testNdop

MCXXFLAGS = \
  -std=c++17 \
  -g \
  -pg \
  -Wno-unused-variable \
  -Wall \
  -Wextra \
  -Werror \
  -pedantic \
  -W \
  -Wshadow \
  -Wdisabled-optimization \
  -march=native \
  -DDEBUG

AM_LDFLAGS = -lggeodesy -lggdatetime

testDorisRinex_SOURCES   = read_rinex.cc
testDorisRinex_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testDorisRinex_LDADD     = $(top_srcdir)/src/libcppdoris.la

testNdop_SOURCES   = test_ndop.cc
testNdop_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testNdop_LDADD     = $(top_srcdir)/src/libcppdoris.la

testAntennaPcv_SOURCES   = test_pcv.cc
testAntennaPcv_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testAntennaPcv_LDADD     = $(top_srcdir)/src/libcppdoris.la
