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

RivetALICEAnalyses.so: $(filter-out $(srcdir)/pluginALICE/tmp*.cc, $(wildcard $(srcdir)/pluginALICE/*.cc))
RivetATLASAnalyses.so: $(filter-out $(srcdir)/pluginATLAS/tmp*.cc, $(wildcard $(srcdir)/pluginATLAS/*.cc))
RivetCDFAnalyses.so: $(filter-out $(srcdir)/pluginCDF/tmp*.cc, $(wildcard $(srcdir)/pluginCDF/*.cc))
RivetCMSAnalyses.so: $(filter-out $(srcdir)/pluginCMS/tmp*.cc, $(wildcard $(srcdir)/pluginCMS/*.cc))
RivetD0Analyses.so: $(filter-out $(srcdir)/pluginD0/tmp*.cc, $(wildcard $(srcdir)/pluginD0/*.cc))
RivetHERAAnalyses.so: $(filter-out $(srcdir)/pluginHERA/tmp*.cc, $(wildcard $(srcdir)/pluginHERA/*.cc))
RivetLEPAnalyses.so: $(filter-out $(srcdir)/pluginLEP/tmp*.cc, $(wildcard $(srcdir)/pluginLEP/*.cc))
RivetLHCbAnalyses.so: $(filter-out $(srcdir)/pluginLHCb/tmp*.cc, $(wildcard $(srcdir)/pluginLHCb/*.cc))
RivetLHCfAnalyses.so: $(filter-out $(srcdir)/pluginLHCf/tmp*.cc, $(wildcard $(srcdir)/pluginLHCf/*.cc))
RivetMCAnalyses.so: $(filter-out $(srcdir)/pluginMC/tmp*.cc, $(wildcard $(srcdir)/pluginMC/*.cc))
RivetMiscAnalyses.so: $(filter-out $(srcdir)/pluginMisc/tmp*.cc, $(wildcard $(srcdir)/pluginMisc/*.cc))
RivetPetraAnalyses.so: $(filter-out $(srcdir)/pluginPetra/tmp*.cc, $(wildcard $(srcdir)/pluginPetra/*.cc))
RivetRHICAnalyses.so: $(filter-out $(srcdir)/pluginRHIC/tmp*.cc, $(wildcard $(srcdir)/pluginRHIC/*.cc))
RivetSPSAnalyses.so: $(filter-out $(srcdir)/pluginSPS/tmp*.cc, $(wildcard $(srcdir)/pluginSPS/*.cc))
RivetTOTEMAnalyses.so: $(filter-out $(srcdir)/pluginTOTEM/tmp*.cc, $(wildcard $(srcdir)/pluginTOTEM/*.cc))

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
