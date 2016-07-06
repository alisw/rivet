SUBDIRS = rivet

if ENABLE_PYEXT

PYEXT_ENV = CC="$(CC)" CXX="$(CXX)" CXXFLAGS="$(PYEXT_CXXFLAGS)" ARCHFLAGS=""

## Always force setup.py, it's not good at guessing what needs to rebuild
all-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --install-lib=build/ --force

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix) --force

uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet
	rm -f  $(DESTDIR)$(RIVET_PYTHONPATH)/rivet-$(VERSION)-py$(PYTHON_VERSION).egg-info

clean-local:
	rm -f $(builddir)/*.pyc
	rm -rf $(builddir)/build
	rm -rf dist

distclean-local:
	rm -f rivet.egg-info

endif
