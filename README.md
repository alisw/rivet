# Welcome

Rivet is a system for preservation of particle-collider analysis logic, analysis
reinterpretation via MC simulations, and the validation and improvement of Monte
Carlo event generator codes. It covers all aspects of collider physics, from
unfolded precision measurements to reconstruction-level searches, and physics
from the Standard Model to BSM theories, and from perturbative jet, boson and
top-quarks to hadron decays, inclusive QCD, and Heavy Ion physics.

Rivet is the most widespread way by which analysis code from the LHC and other
high-energy collider experiments is preserved for comparison to and development
of future theory models. It is used by phenomenologists, MC generator
developers, and experimentalists on the LHC and other facilities. Coding
analyses in Rivet is a great way to publish executable code that extends the
longevity, relevance, and impact of your publications!

These short guides will help you with everything from installation, to first
runs of existing analyses, to writing, running, and plotting results, as well as
how to implement and contribute your own analyses.

Alternatively, feel free to checkout one of the self-guided tutorials (listed towards the
bottom) used for  summer schools and similar. They are suitable for beginners, 
and differently themed depending on the original audience.

Get in touch via the developer mailing list if you need any assistance: [rivet-support@cern.ch](mailto:rivet-support@cern.ch)


## Getting started

[Installation](doc/tutorials/installation.md)

[Rivet via Docker](doc/tutorials/docker.md)

[First rivet run](doc/tutorials/firstrun.md)


## Plotting and run merging

[Plotting with `rivet-mkhtml`](doc/tutorials/plotting.md)

[Customize plots with `make-plots`](doc/tutorials/makeplots.md)

[Merging histograms with `yodamerge` and `rivet-merge`](doc/tutorials/merging.md)


## Advanced running and plotting

[Using analysis options](doc/tutorials/anaoptions.md)

[Preload files, centrality calibration (work in progress)](doc/tutorials/calibration.md)

[Running Rivet on HPC clusters with MPI](doc/tutorials/merging_mpi.md)

[Merging separate physics runs with `rivet-merge` (work in progress)](doc/tutorials/merging2.md)

[Running on a subset of available multiweights](doc/tutorials/multiweights_running.md)

[Fun with multiweights when plotting YODA files (work in progress)](doc/tutorials/multiweights_plotting.md)


## Writing a Rivet analysis
[What is an Analysis?](doc/tutorials/what-analysis.md)

[What is a Projection?](doc/tutorials/projection.md)

[How does Rivet histograms work?](doc/tutorials/rivet-histograms.md)

[Writing a simple analysis](doc/tutorials/simple-analysis.md)

[Writing an analysis with FastJet](doc/tutorials/fastjet.md)

[Contributing a routine](doc/tutorials/anacontrib.md)

[Physics tips and pitfalls](doc/tutorials/tips-pitfalls.md)

[Migrating from Rivet v2 to Rivet v3](doc/tutorials/mig2to3.md)



## Developer topics

[Working with development source (work in progress)](doc/tutorials/developer.md)

[Coding style](doc/tutorials/codingstyle.md)

## Self-guided tutorials

[Resonances, jet physics & weight variations at LHC](doc/tutorials/lhc-basic-tutorial)

[DIS kinematics and final states & analysis options at EIC](doc/tutorials/eic-basic-tutorial)
