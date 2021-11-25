noinst_PROGRAMS = testDorisRinex testAntennaPcv testNdop testDpodExtrapolation testRinexRfoFit testOnePass testTropoH

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

AM_LDFLAGS =-lcppsp3 -lsinex -liers2010 -lggeodesy -lggdatetime

testDorisRinex_SOURCES   = read_rinex.cc
testDorisRinex_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testDorisRinex_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testNdop_SOURCES   = test_ndop.cc
testNdop_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testNdop_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testAntennaPcv_SOURCES   = test_pcv.cc
testAntennaPcv_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testAntennaPcv_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testDpodExtrapolation_SOURCES   = test_pdop_crd.cc
testDpodExtrapolation_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testDpodExtrapolation_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testRinexRfoFit_SOURCES   = test_rnx_fit_rfo.cc
testRinexRfoFit_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testRinexRfoFit_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testOnePass_SOURCES   = test_one_pass.cc
testOnePass_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testOnePass_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)

testTropoH_SOURCES   = test_tropo_h.cpp
testTropoH_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testTropoH_LDADD     = $(top_srcdir)/src/libcppdoris.la $(AM_LDFLAGS)
