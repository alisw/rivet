## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[2.5.1],[rivet@projects.hepforge.org],[Rivet])

## Check and block installation into the src/build dir
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
## Force default prefix to have a path value rather than NONE
if test "$prefix" = "NONE"; then
   prefix=/usr/local
fi

AC_CONFIG_SRCDIR([src/Core/Analysis.cc])
AC_CONFIG_HEADERS([include/Rivet/Config/DummyConfig.hh include/Rivet/Config/RivetConfig.hh include/Rivet/Config/BuildOptions.hh])
AM_INIT_AUTOMAKE([dist-bzip2 -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_CONFIG_MACRO_DIR([m4])
AC_SUBST(LT_OBJDIR)

## Package-specific #defines
AC_DEFINE_UNQUOTED(RIVET_VERSION, "$PACKAGE_VERSION", "Rivet version string")
AC_DEFINE_UNQUOTED(RIVET_NAME, "$PACKAGE_NAME", "Rivet name string")
AC_DEFINE_UNQUOTED(RIVET_STRING, "$PACKAGE_STRING", "Rivet name and version string")
AC_DEFINE_UNQUOTED(RIVET_TARNAME, "$PACKAGE_TARNAME", "Rivet short name string")
AC_DEFINE_UNQUOTED(RIVET_BUGREPORT, "$PACKAGE_BUGREPORT", "Rivet contact email address")

## OS X
AC_CEDAR_OSX

## Work out the LCG platform tag
AC_LCG_TAG

## Set default compiler flags
if test "x$CXXFLAGS" == "x"; then CXXFLAGS="-O2"; fi

## Compiler setup
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

## Store and propagate the compiler identity and flags
RIVETCXX="$CXX"
AC_SUBST(RIVETCXX)
RIVETCXXFLAGS="$CXXFLAGS"
AC_SUBST(RIVETCXXFLAGS)

## Checks for programs.
AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL


## YODA histogramming library
# TODO: we could download, configure, and install YODA automatically... but that is NASTY
AC_CEDAR_LIBRARYANDHEADERS([YODA], , , [AC_MSG_ERROR([YODA is required])])
YODABINPATH=$YODALIBPATH/../bin
AC_SUBST(YODABINPATH)
AC_PATH_PROG(YODACONFIG, yoda-config, [], [$YODALIBPATH/../bin:$PATH])
YODA_PYTHONPATH=""
if test -f "$YODACONFIG"; then
  AC_MSG_CHECKING([YODA version using yoda-config])
  YODA_VERSION=`$YODACONFIG --version`
  AC_MSG_RESULT([$YODA_VERSION])
  YODA_VERSION1=[`echo $YODA_VERSION | cut -d. -f1 | sed -e 's/\([0-9]*\).*/\1/g'`]
  YODA_VERSION2=[`echo $YODA_VERSION | cut -d. -f2 | sed -e 's/\([0-9]*\).*/\1/g'`]
  YODA_VERSION3=[`echo $YODA_VERSION | cut -d. -f3 | sed -e 's/\([0-9]*\).*/\1/g'`]
  let YODA_VERSION_INT=YODA_VERSION1*10000+YODA_VERSION2*100+YODA_VERSION3
  if test $YODA_VERSION_INT -lt 10500; then
    AC_MSG_ERROR([YODA version isn't sufficient: at least version 1.5.0 required])
  fi
  AC_MSG_CHECKING([YODA Python path using yoda-config])
  YODA_PYTHONPATH=`$YODACONFIG --pythonpath`
  AC_MSG_RESULT([$YODA_PYTHONPATH])
fi
AC_SUBST(YODA_PYTHONPATH)


## HepMC event record library
AC_CEDAR_LIBRARYANDHEADERS([HepMC], , , [AC_MSG_ERROR([HepMC is required])])
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$HEPMCINCPATH"
if test -e "$HEPMCINCPATH/HepMC/HepMCDefs.h"; then
  AC_LANG_CONFTEST([AC_LANG_SOURCE([#include <iostream>
#include "HepMC/HepMCDefs.h"
int main() { std::cout << HEPMC_VERSION << std::endl; return 0; }])])
else
  AC_LANG_CONFTEST([AC_LANG_SOURCE([#include <iostream>
#include "HepMC/defs.h"
int main() { std::cout << VERSION << std::endl; return 0; }])])
fi
if test -f conftest.cc; then
  $CXX $CPPFLAGS conftest.cc -o conftest 2>&1 1>&5
elif test -f conftest.C; then
  $CXX $CPPFLAGS conftest.C -o conftest 2>&1 1>&5
else
  $CXX $CPPFLAGS conftest.cpp -o conftest 2>&1 1>&5
fi
hepmc_version=`./conftest`
if test x$hepmc_version != x; then
  let hepmc_major=`echo "$hepmc_version" | cut -d. -f1`
  let hepmc_minor=`echo "$hepmc_version" | cut -d. -f2`
fi
rm -f conftest conftest.cpp conftest.cc conftest.C
HEPMC_VERSION=$hepmc_major$hepmc_minor
AC_MSG_NOTICE([HepMC version is $hepmc_version -> $HEPMC_VERSION])
AC_SUBST(HEPMC_VERSION)
CPPFLAGS=$oldCPPFLAGS


## FastJet clustering library
AC_CEDAR_LIBRARYANDHEADERS([fastjet], , , [AC_MSG_ERROR([FastJet is required])])
AC_PATH_PROG(FJCONFIG, fastjet-config, [], $FASTJETPATH/bin:$PATH)
if test -f "$FJCONFIG"; then
  AC_MSG_CHECKING([FastJet version using fastjet-config])
  fjversion=`$FJCONFIG --version`
  AC_MSG_RESULT([$fjversion])
  fjmajor=$(echo $fjversion | cut -f1 -d.)
  fjminor=$(echo $fjversion | cut -f2 -d.)
  fjmicro=$(echo $fjversion | cut -f3 -d.)
  if test "$fjmajor" -lt 3; then
    AC_MSG_ERROR([FastJet version 3.0.0 or later is required])
  fi
  FASTJETCONFIGLIBADD="$($FJCONFIG --plugins --shared --libs)"
else
  FASTJETCONFIGLIBADD="-L$FASTJETLIBPATH -l$FASTJETLIBNAME"
  FASTJETCONFIGLIBADD="$FASTJETCONFIGLIBADD -lSISConePlugin -lsiscone -lsiscone_spherical"
  FASTJETCONFIGLIBADD="$FASTJETCONFIGLIBADD -lCDFConesPlugin -lD0RunIIConePlugin -lNestedDefsPlugin"
  FASTJETCONFIGLIBADD="$FASTJETCONFIGLIBADD -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin"
  FASTJETCONFIGLIBADD="$FASTJETCONFIGLIBADD -lEECambridgePlugin -lJadePlugin"
fi;
AC_MSG_NOTICE([FastJet LIBADD = $FASTJETCONFIGLIBADD])
AC_SUBST(FASTJETCONFIGLIBADD)
# Check for FastJet headers that require the --enable-all(cxx)plugins option
FASTJET_ERRMSG="Required FastJet plugin headers were not found: did you build FastJet with the --enable-allcxxplugins option?"
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$FASTJETINCPATH"
AC_CHECK_HEADER([fastjet/D0RunIIConePlugin.hh], [], [AC_MSG_ERROR([$FASTJET_ERRMSG])])
AC_CHECK_HEADER([fastjet/TrackJetPlugin.hh], [], [AC_MSG_ERROR([$FASTJET_ERRMSG])])
CPPFLAGS=$oldCPPFLAGS


## GNU Scientific Library
AC_SEARCH_GSL
AC_CEDAR_HEADERS([gsl], , , [AC_MSG_ERROR([GSL (GNU Scientific Library) is required])])
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$GSLINCPATH"
AC_CHECK_HEADER([gsl/gsl_vector.h], [], [AC_MSG_ERROR([GSL vectors not found.])])
CPPFLAGS=$oldCPPFLAGS


dnl ## Boost utility library
dnl BOOST_REQUIRE([1.55.0])
dnl BOOST_FOREACH
dnl BOOST_SMART_PTR
dnl BOOST_FIND_HEADER([boost/lexical_cast.hpp])
dnl BOOST_FIND_HEADER([boost/range.hpp])
dnl BOOST_FIND_HEADER([boost/assign.hpp])
dnl BOOST_BUG_IN_FOREACH


## Disable build/install of standard analyses
AC_ARG_ENABLE([analyses],
  [AC_HELP_STRING(--disable-analyses, [don't try to build or install standard analyses])],
  [], [enable_analyses=yes])
if test x$enable_analyses != xyes; then
   AC_MSG_WARN([Not building standard Rivet analyses, by request])
fi
AM_CONDITIONAL(ENABLE_ANALYSES, [test x$enable_analyses = xyes])


## Enable build/install of unvalidated analyses
AC_ARG_ENABLE([unvalidated],
  [AC_HELP_STRING(--enable-unvalidated, [does nothing, for backward compatibility])],
  [], [enable_unvalidated=yes])
dnl AC_ARG_ENABLE([unvalidated],
dnl   [AC_HELP_STRING(--enable-unvalidated, [build and install unvalidated analyses])],
dnl   [], [enable_unvalidated=yes])
dnl if test x$enable_unvalidated = xyes; then
dnl    AC_MSG_WARN([Building unvalidated Rivet analyses, by request])
dnl else
dnl    AC_MSG_NOTICE([Not building unvalidated Rivet analyses])
dnl fi
dnl AM_CONDITIONAL(ENABLE_UNVALIDATED, [test x$enable_unvalidated = xyes])


## Disable build/install of validated-but-preliminary analyses
AC_ARG_ENABLE([preliminary],
  [AC_HELP_STRING(--disable-preliminary, [does nothing, for backward compatibility])],
  [], [enable_preliminary=yes])
dnl AC_ARG_ENABLE([preliminary],
dnl   [AC_HELP_STRING(--disable-preliminary, [build and install validated-but-preliminary analyses])],
dnl   [], [enable_preliminary=yes])
dnl if test x$enable_preliminary = xyes; then
dnl    AC_MSG_NOTICE([Building preliminary Rivet analyses])
dnl else
dnl    AC_MSG_NOTICE([Not building preliminary Rivet analyses, by request])
dnl fi
dnl AM_CONDITIONAL(ENABLE_PRELIMINARY, [test x$enable_preliminary = xyes])


## Disable build/install of now-obsolete preliminary analyses
AC_ARG_ENABLE([obsolete],
  [AC_HELP_STRING(--disable-obsolete, [does nothing, for backward compatibility])],
  [], [enable_obsolete=yes])
dnl AC_ARG_ENABLE([obsolete],
dnl   [AC_HELP_STRING(--disable-obsolete, [build and install now-obsolete analyses])],
dnl   [], [enable_obsolete=$enable_preliminary])
dnl if test x$enable_obsolete = xyes; then
dnl    AC_MSG_NOTICE([Building obsolete Rivet analyses])
dnl else
dnl    AC_MSG_NOTICE([Not building obsolete Rivet analyses, by request])
dnl fi
dnl AM_CONDITIONAL(ENABLE_OBSOLETE, [test x$enable_obsolete = xyes])


## Build LaTeX docs if possible...
AC_PATH_PROG(PDFLATEX, pdflatex)
AM_CONDITIONAL(WITH_PDFLATEX, [test x$PDFLATEX != x])
## ... unless told otherwise!
AC_ARG_ENABLE([pdfmanual],
  [AC_HELP_STRING(--enable-pdfmanual, [build and install the manual])],
  [], [enable_pdfmanual=no])
if test x$enable_pdfmanual = xyes; then
   AC_MSG_WARN([Building Rivet manual, by request])
fi
AM_CONDITIONAL(ENABLE_PDFMANUAL, [test x$enable_pdfmanual = xyes])

## Build Doxygen documentation if possible
AC_ARG_ENABLE([doxygen],
  [AC_HELP_STRING(--disable-doxygen, [don't try to make Doxygen documentation])],
  [], [enable_doxygen=yes])
if test x$enable_doxygen = xyes; then
   AC_PATH_PROG(DOXYGEN, doxygen)
fi
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])

## Build asciidoc docs if possible
AC_PATH_PROG(ASCIIDOC, asciidoc)
AM_CONDITIONAL(WITH_ASCIIDOC, [test x$ASCIIDOC != x])


## Python extension
AC_ARG_ENABLE(pyext, [AC_HELP_STRING(--disable-pyext,
  [don't build Python module (default=build)])],
  [], [enable_pyext=yes])
AZ_PYTHON_DEFAULT
## Basic Python checks
if test x$enable_pyext == xyes; then
  AZ_PYTHON_PATH
  AZ_PYTHON_VERSION_ENSURE([2.5])
  PYTHON_VERSION=`$PYTHON -c "import sys; print '.'.join(map(str, sys.version_info@<:@:2@:>@));"`
  AC_SUBST(PYTHON_VERSION)
  RIVET_PYTHONPATH=`$PYTHON -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib(prefix='$prefix', plat_specific=True);"`
  AC_SUBST(RIVET_PYTHONPATH)
  ## Test for Python header
  if test -x "$PYTHON"; then
    AC_MSG_CHECKING([for Python include path])
    python_incpath=`$PYTHON -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc();"`
    AC_MSG_RESULT([$python_incpath])
    python_header="$python_incpath/Python.h"
    if test -z "$python_incpath"; then
      AC_MSG_ERROR([Can't build Python extension since include directory cannot be determined from distutils])
      enable_pyext=no
    elif test ! -e "$python_header"; then
      AC_MSG_ERROR([Can't build Python extension since header file $python_header cannot be found])
      enable_pyext=no
    fi
  else
    AC_MSG_ERROR([Can't build Python extension since python can't be found])
    enable_pyext=no
  fi
fi
if test x$enable_pyext == xyes; then
  AC_MSG_NOTICE([All Python build checks successful: 'yoda' Python extension will be built])
fi
AM_CONDITIONAL(ENABLE_PYEXT, [test x$enable_pyext == xyes])


## Set extra Python extension build flags (to cope with Cython output code oddities)
PYEXT_CXXFLAGS="$CXXFLAGS"
AC_CEDAR_CHECKCXXFLAG([-Wno-unused-but-set-variable], [PYEXT_CXXFLAGS="$PYEXT_CXXFLAGS -Wno-unused-but-set-variable"])
AC_CEDAR_CHECKCXXFLAG([-Wno-sign-compare], [PYEXT_CXXFLAGS="$PYEXT_CXXFLAGS -Wno-sign-compare"])
AC_SUBST(PYEXT_CXXFLAGS)


## Cython checks
if test x$enable_pyext == xyes; then
  AM_CHECK_CYTHON([0.18], [:], [:])
  if test x$CYTHON_FOUND = xyes; then
    AC_MSG_NOTICE([Cython >= 0.18 found: Python extension source can be rebuilt (for developers)])
    cython_compiler=$CXX
  fi
fi
AM_CONDITIONAL(WITH_CYTHON, [test x$CYTHON_FOUND = xyes])


if test x$enable_pyext == xyes; then
  AC_MSG_NOTICE([All Python build checks successful: 'rivet' Python extension will be built])
fi
AM_CONDITIONAL(ENABLE_PYEXT, [test x$enable_pyext == xyes])


## Set default build flags
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include"
#AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/include/eigen3"
AM_CPPFLAGS="$AM_CPPFLAGS \$(GSL_CPPFLAGS)"
dnl AM_CPPFLAGS="$AM_CPPFLAGS \$(BOOST_CPPFLAGS)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(YODAINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(HEPMCINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(FASTJETINCPATH)"
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic"])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall"])
AC_CEDAR_CHECKCXXFLAG([-Wno-long-long], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long"])
AC_CEDAR_CHECKCXXFLAG([-Wno-format], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-format"])
dnl AC_CEDAR_CHECKCXXFLAG([-Wno-unused-variable], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-unused-variable"])
AC_CEDAR_CHECKCXXFLAG([-Werror=uninitialized], [AM_CXXFLAGS="$AM_CXXFLAGS -Werror=uninitialized"])
AC_CEDAR_CHECKCXXFLAG([-Werror=delete-non-virtual-dtor], [AM_CXXFLAGS="$AM_CXXFLAGS -Werror=delete-non-virtual-dtor"])


## Debug flag (default=-DNDEBUG, enabled=-g)
AC_ARG_ENABLE([debug], [AC_HELP_STRING(--enable-debug,
  [build with debugging symbols  @<:@default=no@:>@])], [], [enable_debug=no])
if test x$enable_debug == xyes; then
  AM_CXXFLAGS="$AM_CXXFLAGS -g"
fi


## Extra warnings flag (default=none)
AC_ARG_ENABLE([extra-warnings], [AC_HELP_STRING(--enable-extra-warnings,
  [build with extra compiler warnings (recommended for developers)  @<:@default=no@:>@])], [], [enable_extra_warnings=no])
if test x$enable_extra_warnings == xyes; then
   AC_CEDAR_CHECKCXXFLAG([-Wextra], [AM_CXXFLAGS="$AM_CXXFLAGS -Wextra "])
fi


AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_EMPTY_SUBST
AC_CONFIG_FILES(include/Makefile include/Rivet/Makefile)
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(src/Core/Makefile src/Core/yamlcpp/Makefile)
AC_CONFIG_FILES(src/Tools/Makefile)
AC_CONFIG_FILES(src/Projections/Makefile)
AC_CONFIG_FILES(src/Analyses/Makefile)
AC_CONFIG_FILES(test/Makefile)
AC_CONFIG_FILES(pyext/Makefile pyext/rivet/Makefile pyext/setup.py)
AC_CONFIG_FILES(data/Makefile data/refdata/Makefile data/anainfo/Makefile data/plotinfo/Makefile data/texmf/Makefile)
AC_CONFIG_FILES(doc/Makefile)
AC_CONFIG_FILES(doc/rivetversion.sty doc/diffanas)
AC_CONFIG_FILES(bin/Makefile bin/rivet-config bin/rivet-buildplugin)
AC_CONFIG_FILES(Makefile Doxyfile)
AC_CONFIG_FILES(rivetenv.sh rivetenv.csh rivet.pc)

AC_OUTPUT

if test x$enable_pyrivet == xyes; then
   cat <<EOF

************************************************************
RIVET CONFIGURED!

Now build and install (to the $prefix tree) with e.g.
make -j2 && make -j2 install

To use Rivet, we recommend reading HepMC files from a file
or pipe (the latter may be made with mkfifo) using the
'rivet' executable.

For a more pleasant command line experience, you can include
the data/rivet-completion file into your .bashrc file,
or your bash_completion.d directory if you have one.

The rivetenv.*sh files will not be installed, but can help you
to set up a Rivet runtime environment in future.
************************************************************
EOF
fi
