EXTRA_DIST = \
	__init__.py \
    spiresbib.py util.py \
    plotinfo.py aopaths.py \
    core.pyx rivet.pxd \
    core.cpp


if WITH_CYTHON
RIVETINCLUDE = $(top_srcdir)/include/Rivet/
core.cpp: core.pyx rivet.pxd $(RIVETINCLUDE)/Analysis.hh $(RIVETINCLUDE)/AnalysisHandler.hh $(RIVETINCLUDE)/AnalysisLoader.hh $(RIVETINCLUDE)/Run.hh
	cython $(srcdir)/core.pyx --cplus -o $@
else
core.cpp:
	@echo "Not (re)generating core.cpp since Cython is not installed"
endif




## fixes for out-of-source builds, especially "make distcheck"

all-local: fix-out-of-source

FIXSOURCES = $(EXTRA_DIST)

fix-out-of-source: $(FIXSOURCES)
	rm -f $@
	for i in $^; do \
	  orig="$$i"; \
	  build="$$(basename $$i)"; \
	  if [ ! -e $${build} ]; then cp $${orig} $${build}; \
	  echo $${build} >> $@; fi; \
	done
	touch $@

clean-local: fix-out-of-source
	if [ -e $< ]; then if [ -n "$$(cat $<)" ]; then rm -f $$(cat $<); fi; fi
	rm -f $<
