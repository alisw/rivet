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

where the cross-section values are specified in picobarns.
Note that if the HepMC `GenEvent` does not include the cross-section information
and the user also didn't supply a cross-section, the finalised histograms
will have 0 area and the on-the-fly scaling from the above example cannot work.
In this case, it is still possibly to forcibly set a cross-section 
by using a syntax like `ttbar_dilep.yoda.gz:=72.592`. 
To ease scripting and readability,
the syntax `ttbar_dilep.yoda.gz:x72.592` (equivalent to `ttbar_dilep.yoda.gz:72.592`)
for multiplicative scaling will also be accepted.

### Technical note

The final cross-section will be given by the sum of the individual cross-sections.

Note that the `yodamerge` script would really only perform a simple stacking here. 
This is fine when the histogram is normalised to cross-section and the components 
simply just need to be added up, but for a unit-normalised histogram, one would first 
have to undo the divsion by area, add the components with their respective 
cross-section weight, and renormalise the stacked histogram to unity. This is why `yodamerge` 
tends to struggle with files that contain a mix of unit- and cross-section-normalised 
objects, although this can often still be dealt with using a few lines of Python.


## Reentrant safety

In essence, a reentrant-safe analysis is one that, once interrupted mid-run,
has enough information in the output file to be able to pick up where it left off.
The `rivet-merge` script will call the `init()` method of an analysis to book 
its objects into memory, replace all the emoty booked objects with the merged version 
from the combined input files and then run the `finalize()` method once more. 
With that in mind, all fillable objects must be booked in the `init()` method, 
i.e. histograms, profiles, counters. Scatter-type objects can be booked in either 
the `init()` or the `finalize()` method.

The vast majority of analyses shipped with Rivet are reentrant safe. 
Nevertheless, there are a few pitfalls worth highlighting:

 * If an analysis is meant to write out a scatter-type objects, such as
   a ratio or an efficiency, it must book the corresponding numerator 
   and denomintor histograms in the `init()` method in order to be 
   reentrate safe.
 * If the number of histograms booked in the default `init()` method
   depends on the sample (e.g. a beam-energy-dependent booking),
   the routine is not reentrant safe because, for any given run,
   it is not guaranteed that there is always a one-to-one match 
   between the booked analysis objects and the objects written 
   out to file. In such a case, it is preferable to use 
   Rivet's [options mechanism](anaoptions.md) to steer the 
   histogram booking, which ensures that an analysis will
   always be initialised with the right set of objects.
 * It's not possible to use a simple `double` to add up 
   event weights in the main event loop and use the 
   resulting sum to normalise a histogram. There are two reasons: 
   First of all, many Monte Carlo samples will have multiple
   event weights, and so by using a `double` you are probably
   not correctly accounting for the different sum of weights.
   Secondly, when merging files, the main event loop is not 
   executed and so these simple `double` counters will retain
   their initialisation values (usually 0) and the merged result
   cannot receive the intended normalisation. In these cases,
   it is preferable to fill an auxiliary `Counter` object, which
   is then written out to the file, so that it can be loaded back in later.
   The advantage is that this automatically takes care of the multiweights
   and if it's booked with a leading underscore, it will also be ignored
   by the plotting scripts. An example might look like this:
   ```
   CounterPtr mycounter
   book(mycounter, "_auxCounter");
   ...
   normalize(myhisto, *mycounter);
   ```

