SUBDIRS = fjcontrib

noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES = \
  BinnedHistogram.cc \
  Correlators.cc \
  Cuts.cc \
  JetUtils.cc \
  Random.cc \
  Logging.cc \
  ParticleUtils.cc \
  ParticleName.cc \
  Percentile.cc \
  RivetYODA.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  Utils.cc \
  binreloc.c

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DENABLE_BINRELOC -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"

libRivetTools_la_LIBADD = \
  fjcontrib/EnergyCorrelator/libRivetEnergyCorrelator.la \
  fjcontrib/Nsubjettiness/libRivetNsubjettiness.la \
  fjcontrib/RecursiveTools/libRivetRecursiveTools.la
