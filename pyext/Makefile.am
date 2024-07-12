SUBDIRS = rivet

if ENABLE_PYEXT

PYEXT_ENV = CC="$(CXX)" CXX="$(CXX)" CXXFLAGS="$(PYEXT_CXXFLAGS)" ARCHFLAGS="" PYTHONPATH="$(DESTDIR):$(PYTHONPATH)"

all-local:
	cd $(srcdir); \
	$(PYEXT_ENV) RIVET_LOCAL=1 $(PYTHON) $(abs_builddir)/build.py

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(RIVET_PYTHONPATH); \
	cp -r $(abs_builddir)/build/rivet $(DESTDIR)$(RIVET_PYTHONPATH)/

uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet

clean-local:
	rm -rf $(builddir)/build

endif
