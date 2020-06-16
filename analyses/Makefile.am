AUTOMAKE_OPTIONS = -Wno-portability
SHELL := bash -O extglob

if ENABLE_ANALYSES

PLUGIN_DIRS = $(shell ls -d $(srcdir)/plugin*)
PLUGIN_NAMES = $(notdir $(subst plugin,,$(PLUGIN_DIRS)))
PLUGIN_LIBS = $(patsubst %,Rivet%Analyses.so,$(PLUGIN_NAMES))
PLUGIN_INFOFILES = $(shell ls $(abs_srcdir)/plugin*/*.info)
PLUGIN_PLOTFILES = $(shell ls $(abs_srcdir)/plugin*/*.plot)
PLUGIN_YODAFILES = $(shell ls $(abs_srcdir)/plugin*/*.yoda*)
PLUGIN_DATAFILES = $(PLUGIN_INFOFILES) $(PLUGIN_PLOTFILES) $(PLUGIN_YODAFILES)

CLEANFILES = $(PLUGIN_LIBS)
EXTRA_DIST = $(PLUGIN_DIRS)

%.so:
	@+echo && RIVET_BUILDPLUGIN_BEFORE_INSTALL=1 bash $(top_builddir)/bin/rivet-build -j2 $@ $^ -I$(top_builddir)/include

RivetALICEAnalyses.so: $(filter-out $(srcdir)/pluginALICE/tmp*.cc, $(wildcard $(srcdir)/pluginALICE/*.cc))
RivetATLASAnalyses.so: $(filter-out $(srcdir)/pluginATLAS/tmp*.cc, $(wildcard $(srcdir)/pluginATLAS/*.cc))
RivetBABARAnalyses.so: $(filter-out $(srcdir)/pluginBABAR/tmp*.cc, $(wildcard $(srcdir)/pluginBABAR/*.cc))
RivetBELLEAnalyses.so: $(filter-out $(srcdir)/pluginBELLE/tmp*.cc, $(wildcard $(srcdir)/pluginBELLE/*.cc))
RivetBESAnalyses.so: $(filter-out $(srcdir)/pluginBES/tmp*.cc, $(wildcard $(srcdir)/pluginBES/*.cc))
RivetCDFAnalyses.so: $(filter-out $(srcdir)/pluginCDF/tmp*.cc, $(wildcard $(srcdir)/pluginCDF/*.cc))
RivetCESRAnalyses.so: $(filter-out $(srcdir)/pluginCESR/tmp*.cc, $(wildcard $(srcdir)/pluginCESR/*.cc))
RivetCMSAnalyses.so: $(filter-out $(srcdir)/pluginCMS/tmp*.cc, $(wildcard $(srcdir)/pluginCMS/*.cc))
RivetD0Analyses.so: $(filter-out $(srcdir)/pluginD0/tmp*.cc, $(wildcard $(srcdir)/pluginD0/*.cc))
RivetDORISAnalyses.so: $(filter-out $(srcdir)/pluginDORIS/tmp*.cc, $(wildcard $(srcdir)/pluginDORIS/*.cc))
RivetFrascatiAnalyses.so: $(filter-out $(srcdir)/pluginFrascati/tmp*.cc, $(wildcard $(srcdir)/pluginFrascati/*.cc))
RivetHERAAnalyses.so: $(filter-out $(srcdir)/pluginHERA/tmp*.cc, $(wildcard $(srcdir)/pluginHERA/*.cc))
RivetLEPAnalyses.so: $(filter-out $(srcdir)/pluginLEP/tmp*.cc, $(wildcard $(srcdir)/pluginLEP/*.cc))
RivetLHCbAnalyses.so: $(filter-out $(srcdir)/pluginLHCb/tmp*.cc, $(wildcard $(srcdir)/pluginLHCb/*.cc))
RivetLHCfAnalyses.so: $(filter-out $(srcdir)/pluginLHCf/tmp*.cc, $(wildcard $(srcdir)/pluginLHCf/*.cc))
RivetMCAnalyses.so: $(filter-out $(srcdir)/pluginMC/tmp*.cc, $(wildcard $(srcdir)/pluginMC/*.cc))
RivetNovosibirskAnalyses.so: $(filter-out $(srcdir)/pluginNovosibirsk/tmp*.cc, $(wildcard $(srcdir)/pluginNovosibirsk/*.cc))
RivetOrsayAnalyses.so: $(filter-out $(srcdir)/pluginOrsay/tmp*.cc, $(wildcard $(srcdir)/pluginOrsay/*.cc))
RivetMiscAnalyses.so: $(filter-out $(srcdir)/pluginMisc/tmp*.cc, $(wildcard $(srcdir)/pluginMisc/*.cc))
RivetPetraAnalyses.so: $(filter-out $(srcdir)/pluginPetra/tmp*.cc, $(wildcard $(srcdir)/pluginPetra/*.cc))
RivetRHICAnalyses.so: $(filter-out $(srcdir)/pluginRHIC/tmp*.cc, $(wildcard $(srcdir)/pluginRHIC/*.cc))
RivetSLACAnalyses.so: $(filter-out $(srcdir)/pluginSLAC/tmp*.cc, $(wildcard $(srcdir)/pluginSLAC/*.cc))
RivetSPSAnalyses.so: $(filter-out $(srcdir)/pluginSPS/tmp*.cc, $(wildcard $(srcdir)/pluginSPS/*.cc))
RivetTOTEMAnalyses.so: $(filter-out $(srcdir)/pluginTOTEM/tmp*.cc, $(wildcard $(srcdir)/pluginTOTEM/*.cc))
RivetTristanAnalyses.so: $(filter-out $(srcdir)/pluginTristan/tmp*.cc, $(wildcard $(srcdir)/pluginTristan/*.cc))

all-local: $(PLUGIN_LIBS) $(PLUGIN_DATAFILES)
	mkdir -p $(builddir)/data
	for pdir in $(abs_srcdir)/plugin*; do $(LN_S) -f $${pdir}/*.{info,plot,yoda*} $(builddir)/data; done

clean-local:
	rm -rf data

install-exec-local: $(PLUGIN_LIBS)
	$(MKDIR_P) $(DESTDIR)$(libdir)/Rivet
	$(INSTALL) $^ $(DESTDIR)$(libdir)/Rivet

install-data-local: $(PLUGIN_DATAFILES)
	@echo "Installing analysis data files..."
	$(mkdir_p) $(DESTDIR)$(pkgdatadir)
	for pdir in $(builddir)/plugin*; do rsync -aq $${pdir}/*.{info,plot,yoda*} $(DESTDIR)$(pkgdatadir)/ || cp $${pdir}/*.{info,plot,yoda*} $(DESTDIR)$(pkgdatadir)/; done

uninstall-local:
	cd $(DESTDIR)$(libdir) && rm -f $(PLUGIN_LIBS)
	@echo "Uninstalling analysis data files..."
	rm -f $(DESTDIR)$(pkgdatadir)/*.info
	rm -f $(DESTDIR)$(pkgdatadir)/*.plot
	rm -f $(DESTDIR)$(pkgdatadir)/*.yoda*

endif
