# Merging Histograms 101

There are various reasons why one might want to merge together `yoda` files. 
For example, you may want to combine predictions for several sub-processes, 
or outputs from multiple jobs for the same process.
Unfortunately, the devil is in the detail, and simply adding the files
is often not enough to get it right. 
Two utilities are at your disposal to help with the merging:
 * `yodamerge` is a general-purpose script which works for any `yoda`-format file 
 (but has some built-in assumptions), with usage documented below. 
 [\[basic tutorial\]](merging3.md)
 * `rivet-merge` delegates merging of files back to the `Rivet` analyses which produced them.
 [\[basic tutorial\]](merging2.md)
As always, the `--help` flag will also give a lot of information about the respective script
and its limitations. 

## Should I use `yodamerge` or `rivet-merge`?

`yodamerge` is a script built into `yoda` (technically a standalone package from `Rivet`)
which works for any `yoda`-format file. This script does a statistically-correct merging of 
histogram- and profile-type objects. However, when it comes to scatter-type objects,
there are some assumptions/choices which need to be made when merging the `Scatter*D` objects: 
 * should the values of each point simply be added together? 
 (this assumes that each `yoda` file to be merged was generated with the same number of events)
 * should the average be taken for each point? (this assumes that each `yoda` file to be merged was generated with the same number of events)
 * perhaps the points should not be added together, but instead the list of points all `Scatter` objects be concatenated?
 * or finally, one could even just pick the `Scatter` from the first input file and ignore the others.
The answer often depends on the details of the `finalize` method of the parent analysis.
Consider a simple efficiency (a `Scatter2D`) that is constructed from two histograms (`Histo1D` objects).
If only the resulting scatters are written out, the statistical correlations are lost and
it will be impossible to merge the files "correctly". An average might come close, 
but is often not satisfactory.

This is where `rivet-merge` comes in. This script does not make any assumptions about 
how to combine `Scatter*D` objects at all, and instead makes use of re-entrant histogramming,
which starts of combining the pre-finalized versions of the histograms, which are saved to the
output with the prefix `/RAW` prepended to the path. Once combined, the script can then call 
the `finalize()` method of the parent analysis of each analysis object directly in order to correctly 
combine the merged `RAW` histograms into the final `Scatter*D`s.

As a result, `rivet-merge` can **only be used with re-entrant safe routines**. To be re-entrant-safe, 
the `finalize()` method of an analysis should be self-consistent and not book additonal analysis objects:
everything that should end up in the output file must be booked in the initialisation phase.
If you try to merge `yoda` files from non-re-entrant plugins, the script will warn you that the result
will be unpredictable. 

As a rule of thumb, `rivet-merge` is the more sophisticated merging tool, since it has access to the
analysis logic and can actually re-run Rivet over the merged result. 
Please see the [corresponding tutorial](merging2.md) for some examples.
That said, `yodamerge` is good a baseline merging tool that can get far, and in combination with 
a little Python-based post-processing script, anything is possible.
Please see the [corresponding tutorial](merging3.md) for some examples.

