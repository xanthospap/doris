noinst_PROGRAMS = testDorisRinex testSp3Flags testSp3Read testAntennaPcv

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

testSp3Flags_SOURCES   = test_sp3_flags.cc
testSp3Flags_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testSp3Flags_LDADD     = $(top_srcdir)/src/libcppdoris.la

testSp3Read_SOURCES   = test_sp3_read.cc
testSp3Read_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testSp3Read_LDADD     = $(top_srcdir)/src/libcppdoris.la

testAntennaPcv_SOURCES   = test_pcv.cc
testAntennaPcv_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testAntennaPcv_LDADD     = $(top_srcdir)/src/libcppdoris.la
