# Introduction to Rivet

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


## Getting started

[Installation](installation.md)

[Rivet via Docker](docker.md)

[First rivet run](firstrun.md)


## Plotting and run merging

[Plotting with `rivet-mkhtml`](plotting.md)

[Customize plots with `make-plots`](makeplots.md)

[Merging histograms with `yodamerge` and `rivet-merge`](merging.md)


## Advanced running and plotting

[Using analysis options (work in progress)](anaoptions.md)

[Preload files, centrality calibration (work in progress)](calibration.md)

[Merging separate physics runs with `rivet-merge` (work in progress)](merging2.md)

[Running on a subset of available multiweights](multiweights_running.md)

[Fun with multiweights when plotting YODA files (work in progress)](multiweights_plotting.md)


## Writing a Rivet analysis
[What is an Analysis?](what-analysis.md)

[What is a Projection?](projection.md)

[How does Rivet histograms work?](rivet-histograms.md)

[Writing a simple analysis](simple-analysis.md)

[Writing an analysis with FastJet](fastjet.md)

[Contributing a routine](anacontrib.md)

[Migrating from Rivet v2 to Rivet v3](mig2to3.md)


## Developer topics

[Working with development source (work in progress)](developer.md)

[Coding style](codingstyle.md)
