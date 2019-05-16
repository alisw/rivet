dist_noinst_SCRIPTS = \
  mk-analysis-html mk-analysis-latex mk-analysis-json \
  mk-coverage-html mk-coverage-htmls

LATEXSUBSOURCES = \
  preamble.tex intro.tex gettingstarted.tex writinganalyses.tex \
  agilerunmc.tex acknowledgements.tex

DOCSOURCES = $(LATEXSUBSOURCES) \
  compare-histos.txt slip.eps thinker.eps cone.eps bend.eps \
  h-physrev3.bst hepnames.sty JHEP3.cls JHEP.bst make-plots.txt rivet-manual.tex \
  heppennames.sty hepnicenames.sty hepparticles.sty maybemath.sty \
  hepunits.sty underscore.sty microtype.sty

EXTRA_DIST = $(DOCSOURCES) refs.bib

DOCS =
DOCDIRS =


if ENABLE_PYEXT
DOCS += analyses.html
DOCDIRS += analyses
analyses.html: $(top_srcdir)/analyses $(srcdir)/mk-analysis-html
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    RIVET_ANALYSIS_PATH=$(top_builddir)/analyses:: \
    $(PYTHON) $(srcdir)/mk-analysis-html || true
endif


if ENABLE_PYEXT
DOCS += analyses.json
analyses.json: $(top_srcdir)/analyses $(srcdir)/mk-analysis-json
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    RIVET_ANALYSIS_PATH=$(top_builddir)/analyses:: \
    $(PYTHON) $(srcdir)/mk-analysis-json || true
endif


if ENABLE_PYEXT
DOCS += rivet-coverage.html \
  rivet-coverage-heavyiononly.html rivet-coverage-noheavyion.html \
  rivet-coverage-searchesonly.html rivet-coverage-nosearches.html \
  rivet-coverage-nosearches-noheavyion.html
rivet-coverage.html: $(wildcard $(srcdir)/inspire*.json) $(wildcard $(srcdir)/*.rank) $(top_srcdir)/analyses $(srcdir)/mk-coverage-html
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    RIVET_ANALYSIS_PATH=$(top_builddir)/analyses:: \
    $(srcdir)/mk-coverage-htmls || true
endif


if WITH_ASCIIDOC
DOCS += compare-histos.html make-plots.html
EXTRA_DIST += compare-histos.html make-plots.html

compare-histos.html: compare-histos.txt
	asciidoc -a toc compare-histos.txt || true

make-plots.html: make-plots.txt
	asciidoc -a toc make-plots.txt || true
endif


if WITH_PDFLATEX
DOCS += rivet-manual.pdf

if ENABLE_PYEXT
analyses.bib analyses.tex: $(top_srcdir)/analyses $(srcdir)/mk-analysis-latex
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    RIVET_ANALYSIS_PATH=$(top_builddir)/analyses:: \
    $(PYTHON) $(srcdir)/mk-analysis-latex || true
else
analyses.bib analyses.tex:
	> analyses.tex
	> analyses.bib
endif

LATEX = pdflatex --interaction=batchmode
RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"
MAKEIDX = "^[^%]*\\makeindex"

rivet-manual.pdf : rivet-manual.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib
	$(LATEX) $<; true
	egrep -c $(RERUNBIB) rivet-manual.log && (bibtex rivet-manual && cp rivet-manual.toc rivet-manual.toc.bak && $(LATEX) $<); true
	for i in `seq 5`; do if egrep $(RERUN) rivet-manual.log; then echo "LaTeX re-run $i"; cp rivet-manual.toc rivet-manual.toc.bak; $(LATEX) $<; else break; fi; done; true
	if cmp -s rivet-manual.toc rivet-manual.toc.bak; then true; else $(LATEX) $<; true; fi
	rm -f rivet-manual.toc.bak; true

paper.pdf : paper.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib
	$(LATEX) $<; true
	egrep -c $(RERUNBIB) paper.log && (bibtex paper && cp paper.toc paper.toc.bak && $(LATEX) $<); true
	for i in `seq 5`; do if egrep $(RERUN) paper.log; then echo "LaTeX re-run $i"; cp paper.toc paper.toc.bak; $(LATEX) $<; else break; fi; done; true
	if cmp -s paper.toc paper.toc.bak; then true; else $(LATEX) $<; true; fi
	rm -f paper.toc.bak; true

endif


################


RSH = rsync
DEST = login.hepforge.org:rivet/public_html/


.PHONY = doc upload arxivtar

doc: $(DOCS)
	@true

upload: $(DOCS) $(DOCDIRS)
	$(RSH) -r $? $(DEST)

arxivtar: $(DOCSOURCES) rivet-manual.bbl
	for i in *.png; do convert $$i $${i/.png/.eps}; done
	tar czf rivet-manual.tar.gz \
rivet-manual.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib rivetversion.sty rivet-manual.bbl \
hepnames.sty hepnicenames.sty hepparticles.sty heppennames.sty hepunits.sty maybemath.sty microtype.sty underscore.sty \
bend.eps cone.eps thinker.eps slip.eps \
h-physrev3.bst JHEP3.cls JHEP.bst


all-local:
	@echo "No default target in doc: run 'make doc' to build documentation"
	@true

install-data-local:
	@echo "No local docs installation: copy the documentation to wherever you require"
	@true


mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf rivet-coverage-*.html $(DOCS) $(DOCDIRS)
