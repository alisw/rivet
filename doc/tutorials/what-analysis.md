## What is a "Rivet analysis"?
Before starting to write your own Rivet analysis, it might be beneficial to read a little more about
what an analysis is on the technical side. Briefly:

There is one Analysis for every physics paper implemented.

They produce the histograms which can be compared to the published plots in the paper.

They calculate event properties and implement kinematic cuts using [Projections](projections.md).

They book, fill and output histograms using the YODA interface via auxilliary RivetHistogramming code. Each event is passed to the analysis and operated on by the projections. The result of the projections operating on the event determines whether the event is accepted and what is plotted.

All the Analyses to be used are instantiated at the beginning of a Rivet run.

All Analyses inherit from the Analysis.hh abstract base class.

### Writing your own analysis
Rivet uses ''pluggable'' analyses: this system allows users to build and run their own analysis without needing to re-build the Rivet library.

See [Writing a simple analysis](simple-analysis.md) for instructions on writing your own analysis routine.

### How analyses get loaded

When you build your new analysis a `RivetMyAnalysis.so` library will be created, which Rivet needs to find.
Rivet scans several directories at runtime to find the library files containing analyses. Specifically, the following search locations are tried, in order:

 * a directory listed in `RIVET_ANALYSIS_PATH}
 * the directory where `libRivet.so` is installed (i.e. the Rivet `$prefix/lib` directory)
 * the current directory.

If a duplicate analysis is found in more than one location, Rivet will complain but will not crash. The first version to be found will be used.

Note that (to reduce the number of attempted loadings) *the library name must contain the word `Rivet` and end in the appropriate shared library suffix for the OS: this is .so for Linux and Macs*.

### Real Data

YODA files for all measurements are distributed with Rivet for comparison to the generated plots. These files are obtained from [HepData](http://www.hepdata.net), but distributed with the Rivet code for standalone running.

These data files are also used to auto-book the histograms which are filled by the Monte Carlo. See [Book and use histograms](rivet-histograms.md).
