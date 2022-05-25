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

if ENABLE_HEPMC_3
libRivetTools_la_SOURCES += RivetHepMC_3.cc ReaderCompressedAscii.cc WriterCompressedAscii.cc
else
libRivetTools_la_SOURCES += RivetHepMC_2.cc
endif

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) -DENABLE_BINRELOC -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"
