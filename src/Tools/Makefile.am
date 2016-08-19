noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES = \
  BinnedHistogram.cc \
  Cuts.cc \
  Logging.cc \
  ParticleName.cc \
  RivetYODA.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  Utils.cc \
  binreloc.c \
  mt2_bisect.cc

dist_noinst_HEADERS = binreloc.h mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DENABLE_BINRELOC -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"
