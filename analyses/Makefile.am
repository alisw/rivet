AUTOMAKE_OPTIONS = -Wno-portability
SHELL := bash -O extglob

if ENABLE_ANALYSES

PLUGIN_DIRS = $(shell ls -d $(srcdir)/plugin*)
PLUGIN_NAMES = $(notdir $(subst plugin,,$(PLUGIN_DIRS)))
PLUGIN_LIBS = $(patsubst %,Rivet%Analyses.so,$(PLUGIN_NAMES))
PLUGIN_DATAFILES = $(shell ls $(abs_srcdir)/plugin*/*.{info,plot,yoda})

CLEANFILES = $(PLUGIN_LIBS)
EXTRA_DIST = $(PLUGIN_DIRS)

%.so:
	@+echo && RIVET_BUILDPLUGIN_BEFORE_INSTALL=1 bash $(top_builddir)/bin/rivet-buildplugin -j2 $@ $^ -I$(top_builddir)/include

RivetALICEAnalyses.so: $(shell ls $(srcdir)/pluginALICE/!(tmp)*.cc)
RivetATLASAnalyses.so: $(shell ls $(srcdir)/pluginATLAS/!(tmp)*.cc)
RivetCDFAnalyses.so  : $(shell ls $(srcdir)/pluginCDF/!(tmp)*.cc)
RivetCMSAnalyses.so  : $(shell ls $(srcdir)/pluginCMS/!(tmp)*.cc)
RivetD0Analyses.so   : $(shell ls $(srcdir)/pluginD0/!(tmp)*.cc)
RivetHERAAnalyses.so : $(shell ls $(srcdir)/pluginHERA/!(tmp)*.cc)
RivetLEPAnalyses.so  : $(shell ls $(srcdir)/pluginLEP/!(tmp)*.cc)
RivetLHCbAnalyses.so : $(shell ls $(srcdir)/pluginLHCb/!(tmp)*.cc)
RivetLHCfAnalyses.so : $(shell ls $(srcdir)/pluginLHCf/!(tmp)*.cc)
RivetMCAnalyses.so   : $(shell ls $(srcdir)/pluginMC/!(tmp)*.cc)
RivetMiscAnalyses.so : $(shell ls $(srcdir)/pluginMisc/!(tmp)*.cc)
RivetPetraAnalyses.so: $(shell ls $(srcdir)/pluginPetra/!(tmp)*.cc)
RivetRHICAnalyses.so : $(shell ls $(srcdir)/pluginRHIC/!(tmp)*.cc)
RivetSPSAnalyses.so  : $(shell ls $(srcdir)/pluginSPS/!(tmp)*.cc)
RivetTOTEMAnalyses.so: $(shell ls $(srcdir)/pluginTOTEM/!(tmp)*.cc)

all-local: $(PLUGIN_LIBS) $(PLUGIN_DATAFILES)
	mkdir -p $(builddir)/data
	$(LN_S) -f $(abs_srcdir)/plugin*/*.{info,plot,yoda} $(builddir)/data

clean-local:
	rm -rf data

install-exec-local: $(PLUGIN_LIBS)
	$(INSTALL) $^ $(DESTDIR)$(libdir)

install-data-local: $(PLUGIN_DATAFILES)
	@echo "Installing analysis data files..."
	$(mkdir_p) $(DESTDIR)$(pkgdatadir)
	cp $(abs_srcdir)/plugin*/*.{info,plot,yoda} $(DESTDIR)$(pkgdatadir)/
#@$(install_sh_DATA) $^ $(DESTDIR)$(pkgdatadir)

uninstall-local:
	cd $(DESTDIR)$(libdir) && rm -f $(PLUGIN_LIBS)
	@echo "Uninstalling analysis data files..."
	rm -f $(DESTDIR)$(pkgdatadir)/*.{info,plot,yoda}

endif
