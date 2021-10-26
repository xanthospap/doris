noinst_PROGRAMS = testDorisRinex testAntennaPcv testNdop testDpodExtrapolation testRinexRfoFit

MCXXFLAGS = \
  -std=c++20 \
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

AM_LDFLAGS = -lsinex -lggeodesy -lggdatetime

testDorisRinex_SOURCES   = read_rinex.cc
testDorisRinex_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testDorisRinex_LDADD     = $(top_srcdir)/src/libcppdoris.la

testNdop_SOURCES   = test_ndop.cc
testNdop_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testNdop_LDADD     = $(top_srcdir)/src/libcppdoris.la

testAntennaPcv_SOURCES   = test_pcv.cc
testAntennaPcv_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testAntennaPcv_LDADD     = $(top_srcdir)/src/libcppdoris.la

testDpodExtrapolation_SOURCES   = test_pdop_crd.cc
testDpodExtrapolation_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testDpodExtrapolation_LDADD     = $(top_srcdir)/src/libcppdoris.la

testRinexRfoFit_SOURCES   = test_rnx_fit_rfo.cc
testRinexRfoFit_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testRinexRfoFit_LDADD     = $(top_srcdir)/src/libcppdoris.la