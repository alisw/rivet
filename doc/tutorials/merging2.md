# Basic `rivet-merge` tutorial

As of version 2.7 the concept of a *reentrant analysis* was introduced in Rivet.
This means that an analysis is written so that the `finalize` method can be run
several times. This also means that `rivet` can read in a previously produced
`yoda` file and re-run `finalize`. One of the advantages of this is that it is
then possible to merge several runs into one in a statistically correct way
using the `rivet-merge` script.

In the following, we will distinguish _equivalent merging_,
where all output files are assumed to come from a CPU-parallelised run
over the same pool of generated events using different random seeds, 
and _non-equivelent merging_ (or _stacking_), where the files are assumed to 
correspond to different processes and are just meant to be stacked.


## Equivalent merging

Assuming you have produced files named `run00.yoda`, ...,
`run99.yoda` with different random seeds for the event generator, you simply do
```
rivet-merge -e run[0-9][0-9].yoda -o sumrun.yoda
```
and `sumrun.yoda` will contain the sum of all runs. (Note the `-e` flag that
indicates that the given `yoda` files are equivalent.)

### Technical note

In most cases, a histogram will be scaled by a factor in the `finalize` method,
either cross-section / sum of weights or 1 / area. In practice, this scaling
needs to be "undone" for every input histogram first. The unscaled histograms 
are then stacked and the individual scaling factors combined by summing the
inverse input scale factors and taking the inverse of the sum.
The combined scaling factor is then applied to the stacked histogram.

For the final cross-section, a weighted average will be constructed from
the per-file cross-sections using the respective sum of weights in a given
file as the weight.


## Stacking

Let's assume that the analysis measures a jet spectrum over a large range 
of transverse momenta. It is then difficult to get reasonable statistics 
for all transverse momenta in a single generator run. However, the event generator 
most likely have the possibility to select upper and lower cuts on the 
transverse momentum of the hard partonic sub-process. Although such cuts
do not correspond to clean cut in the jet transverse momenta, we can still merge
them together with `rivet-merge`. If you have four files produced with cuts on
the hard sub-process transverse momenta as indicated by the file names
`run_10-20.yoda`, `run_20-50.yoda`, `run_50-100.yoda`, and `run_100-inf.yoda`,
they are simply merged together with
```
rivet-merge run_10-20.yoda run_20-50.yoda run_50-100.yoda run_100-inf.yoda 
```
where the merged result is found in the default
output file `Rivet.yoda`. Note the absence of the `-e` flag.


Apart from slicing, one could also envisage that a process is split based
on final state. For instance, consider top-quark pair production in the dileptonic, 
semi-leponic and all-hadronic channels, with corrsponding `yoda` files 
`ttbar_dilep.yoda.gz`, `ttbar_singlep.yoda.gz` and `ttbar_hadronic.yoda.gz` respectively.
Let's assume the generator-level cross-section was not known when the files
were produced and so the cross-section was manually set to unity at the time.
It is possible to stack the files and pass the cross-section on the fly as follows:
```
rivet-merge -o ttbar.yoda.gz ttbar_dilep.yoda.gz:72.592 ttbar_singlep.yoda.gz:302.06 ttbar_hadronic.yoda.gz:314.12
```

### Technical note

The final cross-section will be given by the sum of the individual cross-sections.

Note that the `yodamerge` script would really only perform a simple stacking here. 
This is fine when the histogram is normalised to cross-section and the components 
really just need to be added up, but for a unit-normalised histogram, one would first 
have to undo the divsion by area, add the components with their respective 
cross-section weight, and renormalise the stack to unity. This is why `yodamerge` 
tends to struggle with files that contain a mix of unit- and cross-section-normalised 
objects, although this can often still be dealt with using a few lines of Python.

