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
  wget https://gitlab.com/hepcedar/rivetbootstrap/raw/3.1.0/rivet-bootstrap
  chmod +x rivet-bootstrap
```
(Replace the version string as appropriate if you want to install other versions of Rivet.)

2. **Check the options.** Look at the header of the script to see all variables which you can set, e.g. to skip installation of certain dependencies if they are available in your system:
```
  less rivet-bootstrap ## and read...
```


3. **Run the script.** By default it will install to `$PWD/local`, where `$PWD` is the current directory. If you need to change that, specify the corresponding values on the command line. Examples:
```
./rivet-bootstrap
  # or, e.g.
INSTALL_PREFIX=$HOME/software/rivet MAKE="make -j8" ./rivet-bootstrap
```
We will refer to the installation root path as `$PREFIX`.

*Other variables used in the script, such as version numbers, can be overridden
 this way. Those other than `INSTALL_PREFIX` and `MAKE` are mainly of interest
 to developers, though. It is possible that you will need a developer-style
 setup of the Rivet development version rather than a released version: in this
 case you will need to have essential developer tools like autotools and Cython
 already installed, then pass `INSTALL_RIVETDEV=1` as a command-line variable.*


## Setting up the environment

After the script grinds away for a while, it will tell you that it is finished and how to set up a runtime environment (similar to that used inside the build script) for running Rivet. A sourceable rivetenv.(c)sh script is provided for (c)sh shell users to help set up this environment. Here's how to set up the environment and then test the `rivet` program's help feature and analysis listing:

```
  source $PREFIX/rivetenv.sh
  rivet --help
  rivet --list-analyses
```

If that works, everything is installed correctly. If you are using the `bash` shell in your terminal, then Rivet will offer you programmable tab completion: try typing {{{rivet}}} and pressing the Tab key!

**You may wish to add the environment variable settings to your `~/.bashrc` shell config file, so that Rivet will work without needing any special session setup.**


## Troubleshooting

The main things that can go wrong during Rivet installation and setup are
related to the Python interface. Rivet relies heavily on Python interfaces to
the C++ libraries for its user-facing command-line scripts, but this requires
a suitable Python environment before installation.

If you hit problems, a good first idea is to make sure that the installation
prefix is treated as a valid location for Python to load packages from. If
installing into a system location like `/usr/local` this should be automatic; if
doing something more custom, you may need to set or append to the `PYTHONPATH`
environment variable. The usual Python install path is
`$INSTALL_PREFIX/lib/pythonX.Y/site-packages`, but it's best to check your
system for variations; once confirmed, put
`export PYTHONPATH=$INSTALL_PREFIX/lib/pythonX.Y/site-packages:$PYTHONPATH`
into your shell environment, perhaps permanently into your `.bashrc` file.

The second issue is specific to Macs, on which the system Python is configured
to refuse to load modules that load "untrusted" binary libraries... like `libRivet`.
We have not found a satisfactory workaround for this restriction, which seems
designed under the assumption that Mac users will not be developers of Python
extensions to compiled non-system libraries. Since this is unrealistic for
HEP use, we recommend avoiding the Mac system Python and compiler suite entirely,
in favour of developer packages from e.g. Homebrew or Conda.
