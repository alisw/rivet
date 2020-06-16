# Merging histograms with `yodamerge` and `rivet-merge` 

There are various reasons why one might want to merge together `yoda` files. For example, you may want to combine predictions for several sub-processes, or outputs from multiple jobs for the same process.
Two utilities are at your disposal to help do this merging: `yodamerge` and `rivet-merge`.
`yodamerge` is a general-purpose script which works for any `yoda`-format file (but has some built-in assumptions), with usage documented below.
`rivet-merge` delegates merging of files back to the `Rivet` analyses which produced them, and is documented in more detail [here](merging2.md).
For a discussion on the differences between the two scripts, check the [dedicated section](https://gitlab.com/hepcedar/rivet/-/blob/doc-merging/doc/tutorials/merging.md#should-i-use-yodamerge-or-rivet-merge). 

First of all, to get information about `yodamerge`, try to run:
```
yodamerge --help
```
which is a treasure-trove of information. The below is just a basic tutorial with simple examples.

# Basic examples of using `yodamerge`

Note that example files `file1.yoda` `file2.yoda` are included in this directory.

The basic syntax of `yodamerge` is:
```
yodamerge -o newyoda.yoda file1.yoda file2.yoda # file3.yoda... etc
```
This creates a new file `newyoda.yoda` in which any `Histo*D`, `Profile*D`, and `Counter` objects are correctly statistically combined, by undoing any scaling which has been applied, then summing the objects and re-apply the required scaling. 
```
yodamerge --add -o newyoda.yoda file1.yoda file2.yoda # file3.yoda... etc
```
This does the same as the previous example, but forces simple stacking of `Histo` objects and other types, without applying the unscaling-rescaling procedure described (briefly) above. This will literally just add the numbers in each bin of `Histo` objects, and not worry that one run may be larger than the other (normally a weighted average is taken). .

```
yodamerge --add -o newyoda.yoda file1.yoda:1.23 file2.yoda:4.56 # file3.yoda:7.89... etc
```
Another useful feature shown above is that one can specify a custom scaling of the individual files, using the `<file name>:<scaling float>` convention.

In all these cases, assumption have been made here for `Scatter*D` objects. See below for details of how to modify those assumptions.

# Options for merging `Scatter*D` objects

You may have noticed the warning:
```
WARNING: Scatter2D /REF/ATLAS_2019_I1725190/d01-x01-y01 merge assumes asymptotic statistics and equal run sizes
WARNING: Scatter2D /REF/ATLAS_2019_I1725190/d02-x01-y01 merge assumes asymptotic statistics and equal run sizes
```
Which is printed for every `Scatter*D` object in the files.

This is because the default behaviour for merging scatters is that the average of each point of a `Scatter` object is taken. This seemingly straightforward behaviour actually assumes, technically, that the `yoda` files which are being merged had the same number of events processed which they were filled. This default behaviour is equivalent to using the `assume_mean` argument for the `--s2d-mode` option (for 2D `Scatter` objects). For this and all that follows, similar options exist for 1D and 3D objects (`--s1d-mode` and `--s3d-mode`).
```
yodamerge file1.yoda file2.yoda -o newyoda.yoda --s2d-mode assume_mean
```
Other choices are possible... one may want to simply add the points from the `Scatter*D` objects together without averaging. This can be achieved using `--s2d-mode add`:

```
yodamerge file1.yoda file2.yoda -o newyoda.yoda --s2d-mode add 
```

Another behaviour might be that the points should not be added or averaged at all, but instead that the list of points from the various `Scatter*D` objects should be concatenated. So two 6-point `Scatters` would be merged into a single 12-point `Scatter`. This would be achieved using the `combine` argument :

```
yodamerge file1.yoda file2.yoda -o newyoda.yoda --s2d-mode combine 
```

Or finally, it may simply be that it does not make much sense to merge `Scatter` objects and that the safest thing to do is pick one of them, the fist to be encountered, and write that out to the final `yoda` file. In this case, use `--s2d-mode first`

# Other options

There are a few other options which may be of use:

```
--type-mismatch-mode MODE # 'first'or 'scatter'
```
This option may be of use if you encounter objects with the same name but different types. Either you can just pick the first one and skip the rest, or convert everything to a `Scatter` object and treat it that way.


```
--no-veto-empty       
```
By default, empty (sumW=0) data objects are not written out to the final `yoda` file. This option writes out such objects.

# Do it yourself!

If you need features which are not implemented here, or have your own specific mergign requirements, then you can make a local copy of `yodamerge` and use the very powerful `Python` API to write in your own merging prescriptions. 

## Should I use `yodamerge` or `rivet-merge`?

`yodamerge` is a script built into `yoda` (which is technically a standalone package from `Rivet`) which works for any `yoda`-format file. This script does a statistically-correct merhing of Histogram and Profile objects... But there are some assumptions/choices which need to be made when merging `Scatter*D` objects: 
 * should the values of each point simply be added together ? (this assumes that each `yoda` file to be merged was generated with the same number of events)
 * should the average be taken for each point ? (this assumes that each `yoda` file to be merged was generated with the same number of events)
 * perhaps the points should not be added together, but instead the list of points all `Scatter` objects be concatenated?
 * or finally, one could even just pick the `Scatter` from the first input file and ignore the others.

Knowing how to combine objects like `Scatter*D`s correctly may be analysis-specific, and may also require one to have access to the history of how the `yoda` file was filled. This is where `rivet-merge` comes in. This script does not make any assumptions at all about how to combine `Scatter*D` objects, and instead uses to re-entrant histogramming feature of `Rivet` (which keeps the 'history' of how analysis objects were filled in the form of `RAW` histograms). It can then call directly the `finalize()` method of the parent analysis of each analysis object to correctly combine the merged `RAW` histograms into the final `Scatter*D`s. One final assumption needs to be made, which is whether the runs you are merging are for identical processes (in which case use the `-e` or `--equiv` option) or not. 

The corrollary of the fact that `rivet-merge` uses re-entrant-safe `Rivet`-analysis `finalize()` methods is that `rivet-merge` can **only** be used for `yoda` files which result from re-entrant-safe `Rivet` analyses. If you try to merge `yoda` files from non-re-entrant plugins, the script will warn you, and then you should use `yodamerge` instead.

`rivet-merge` is documented on a separate page: [Merging separate physics runs with `rivet-merge`](merging2.md)
