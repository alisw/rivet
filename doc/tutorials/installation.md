# Installation

The easiest way to start using Rivet is via the Docker system (see https://www.docker.com/), which is like a lightweight Linux virtual machine that's easily installed and run on Mac, Linux and Windows machines. See the dedicated [Docker](docker.md) tutorial page for more detail.

The rest of these instructions are mainly aimed at users who want to natively install and run a release of Rivet on their machine.


## Native installation of Rivet and all dependencies

The simplest Rivet installation uses a "bootstrap" script to install Rivet and all its dependencies from release tarballs.


### Prerequisites

Python header files are required. On Ubuntu and other Debian Linux derivatives,
you can use this command to install the necessary files system-wide: `sudo apt-get install python-dev`.
You will also need a C++ compiler capable of building the C++14 dialect: on systems like
CERN lxplus you can get such an environment (and also a fixed LaTeX for plotting) with:
```
source /cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc62-opt/Python-env.sh
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2016/bin/x86_64-linux:$PATH
```

### Installation

1. **Download the bootstrap script** into a temporary working directory, and make it executable:
```
  cd /scratch/rivet
  wget https://gitlab.com/hepcedar/rivetbootstrap/raw/3.1.3/rivet-bootstrap
  chmod +x rivet-bootstrap
```
(Replace the version string as appropriate if you want to install other versions of Rivet.)

2. **Check the options.** Look at the header of the script to see all variables which you can set, e.g. to skip installation of certain dependencies if they are available in your system:
```
  less rivet-bootstrap ## and read...
```


3. **Run the script.** By default it will install the whole suite of Rivet dependencies
and Rivet itself to `$PWD/local`, where `$PWD` is the current directory.
We will refer to this installation root path as `ROOT`: where the word `ROOT`
appears below, in commands for you to type in, you should replace it with the actual
installation prefix you used when running `rivet-bootstrap`

If you
need to change that, specify the corresponding values on the command line. Other variables
used by the script can be set at the same time if you wish. Examples:
```
# To install to $PWD/local:
./rivet-bootstrap
```
or
```
# To install to ~/software/rivet
INSTALL_PREFIX=$HOME/software/rivet MAKE="make -j8" ./rivet-bootstrap
```

*Other variables used in the script, such as version numbers, can be overridden
 this way. Those other than `INSTALL_PREFIX` and `MAKE` are mainly of interest
 to developers, though. It is possible that you will need a developer-style
 setup of the Rivet development version rather than a released version: in this
 case you will need to have essential developer tools like autotools and Cython
 already installed, then pass `INSTALL_RIVETDEV=1` as a command-line variable.*


## Setting up the environment

After the script grinds away for a while, it will tell you that it is finished and how to set up a runtime environment (similar to that used inside the build script) for running Rivet. A sourceable rivetenv.(c)sh script is provided for (c)sh shell users to help set up this environment. Here's how to set up the environment and then test the `rivet` program's help feature and analysis listing:

```
  source ROOT/rivetenv.sh
  rivet --help
  rivet --list-analyses
```

If that works, everything is installed correctly. If you are using the `bash` shell in your terminal, then Rivet will offer you programmable tab completion: try typing `rivet ` and pressing the `Tab` key! Note the space before pressing `Tab`: if you just type `rivet` + `Tab`, without a space, you'll get a handy list of all the Rivet tool programs.

**You may wish to add the environment variable settings to your `~/.bashrc` shell config file, so that Rivet will work without needing any special session setup.**


## Troubleshooting

### PYTHONPATH and other path-search issues

The main things that can go wrong during Rivet installation and setup are
related to the Python interface. Rivet relies heavily on Python interfaces to
the C++ libraries for its user-facing command-line scripts, but this requires
a suitable Python environment before installation.

If you hit problems, a good first idea is to make sure that the installation
prefix is treated as a valid location for Python to load packages from *before*
running the installation process. If installing into a system location like
`/usr/local` this should be automatic, but if doing something more custom, you
may need to set or append to the `PYTHONPATH` environment variable. (This, and
equivalent additions to the `PATH` and `(DY)LD_LIBRARY_PATH` environment variables
for executable and shared-library loading respectively, is done post-build
by the `rivetenv.(c)sh` setup script.)

The usual install path for Rivet's Python package is in
`$INSTALL_PREFIX/lib/pythonX.Y/site-packages`, but it's best to check your
system for variations. Once confirmed, put
`export PYTHONPATH=$INSTALL_PREFIX/lib/pythonX.Y/site-packages:$PYTHONPATH`
into your shell environment, and then run the installer.

### Mac untrusted binaries

The second issue is specific to Macs, on which the system Python is configured
to refuse to load modules that load "untrusted" binary libraries... like `libRivet`.
We have not found a satisfactory workaround for this restriction, which seems
designed under the assumption that Mac users will not be developers of Python
extensions to compiled non-system libraries. Since this is unrealistic for
HEP use, we recommend avoiding the Mac system Python and compiler suite entirely,
in favour of developer packages from e.g. Homebrew or Conda.
