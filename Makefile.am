ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext data include bin analyses test doc

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-yoda=$(YODAPATH) --with-hepmc=$(HEPMCPATH) --with-fastjet=$(FASTJETPATH)

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = rivet.pc

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox pyclean

clean-local:
	@rm -rf a.out

install-data-local:
	cd doc && $(MAKE) dat json
	$(mkdir_p) $(DESTDIR)$(pkgdatadir)
	$(install_sh_DATA) doc/analists/analyses.dat $(DESTDIR)$(pkgdatadir)/
	$(install_sh_DATA) doc/analists/analyses.json $(DESTDIR)$(pkgdatadir)/

## Upload tarballs to HepForge
RSH = rsync
DEST = login.hepforge.org:rivet/downloads/
upload: dist
	$(RSH) $(DIST_ARCHIVES) $(DEST)
