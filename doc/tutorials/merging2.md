# Merging rivet runs with `rivet-merge`

As of version 2.7 the concept of a *reentrant* analysis was introduced in Rivet.
This means that an analysis is written so that the `finalize` method can be run
several times. This also means that `rivet` can read in a previously produced
`yoda` file and re-run `finalize`. One of the advantages of this is that it is
then possible to merge several runs into one in a statistically correct way
using the `rivet-merge` script.

The `rivet-merge` script addresses two different scenaria. One casenis when you
eg. have made the same generator run spread our over many CPUs (with different
random seeds) obtaining several equivalent `yoda` files. The other is when you
eg. have run the same analysis for two different processes and want to merge
them together.

The first case is easiest, and can in in most cases be handled just as well with
`yodamerge` script. Assuming you have produced files named `run00.yoda`, ...,
`run99.yoda` woth different random seeds for the event generator, you simply do
```
rivet-merge -e run[0-9][0-9].yoda -o sumrun.yoda
```
and `sumrun.yoda` will contain the sum of all runs. (note the `-e` flag that
indicates that the given `yoda` files are equivalent)

For the second case we assume that the analysis measures a jet spectrum over
a large range of transverse momenta. It is then difficult to get reasonable
statistics for all transverse momenta in a single generator run. However, the
event generator most likely have the possibility to select upper and lower cuts
on the transverse momentum of the hard partonic sub-process. Although such cuts
do not correspond to clean cut in the jet transverse momenta, we can still merge
them together with `rivet-merge`. If you have four files produced with cuts on
the hard sub-process transverse momenta as indicated by the file names
`run_10-20.yoda`, `run_20-50.yoda`, `run_50-100.yoda`, and `run_100-inf.yoda`,
they are simply merged together with
```
rivet-merge run_10-20.yoda run_20-50.yoda run_50-100.yoda run_100-inf.yoda 
```
(note the absence of `-e`) where the merged result is found in the default
output file `Rivet.yoda`.

