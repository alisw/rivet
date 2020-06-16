# How to contribute a routine

**We are in the process of streamlining the routine contribution procedure
using a merge-request template as well as a dedicated CI pipeline.
In the meantime, feel free to make a merge request (if you know where
the relevant files should go) or just open an issue on our issue tracker.**


1.  Make sure the reference data file is fully in sync with what's on HEPData.

2.  Please provide the routine `*.cc`, `*.plot` and `*.info` files. These will
    go into the relevant `analyses/pluginEXPERIMENT`directory.

3.  Please provide the `yoda` file used to produce the validation plots.
    Note that it's not necessary to actually tar up the `rivet-mkhtml` 
    output, the yoda file used to produce this is fully sufficient.

4.  Please identify a suitable sample in our 
    [suite of HepMC validation files](https://rivetval.web.cern.ch/rivetval/HEPMC/)
    and specify its name in the info file e.g. like so
    ```
    ReleaseTests:
     - $A LHC-13-Top-L
    ```
    These files only have 1000 events and are mainly used for numerical regression tests.
    
