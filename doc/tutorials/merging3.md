# Basic `yodamerge` tutorial

Note that example files `file1.yoda` `file2.yoda` are included in this directory.

The basic syntax of `yodamerge` is:
```
yodamerge -o newyoda.yoda file1.yoda file2.yoda # file3.yoda... etc
```
This creates a new file `newyoda.yoda` in which any `Histo*D`, `Profile*D`, and `Counter` objects are correctly statistically
combined, by undoing any scaling which has been applied, then summing the objects and re-applying the required scaling. 

```
yodastack -o newyoda.yoda file1.yoda file2.yoda # file3.yoda... etc
```
This does the same as the previous example, but forces simple stacking of histogram objects and other types, 
without applying the unscaling-rescaling procedure described (briefly) above. This will literally just add 
the numbers in each bin of histogram objects, 
and not worry that one run may be larger than the other (normally a weighted average is taken).

```
yodastack -o newyoda.yoda file1.yoda:1.23 file2.yoda:4.56 # file3.yoda:7.89... etc
```
Another useful feature shown above is that one can specify a custom scaling of the individual files, 
using the `<file name>:<scaling float>` convention.

In all these cases, assumption have been made here for `Scatter*D` objects. See below for details of how to modify those assumptions.

## Options for merging `Scatter*D` objects

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

## Other options

There are a few other options which may be of use:

```
--type-mismatch-mode MODE # 'first'or 'scatter'
```
This option may be of use if you encounter objects with the same name but different types. Either you can just pick the first one and skip the rest, or convert everything to a `Scatter` object and treat it that way.


```
--no-veto-empty       
```
By default, empty (sumW=0) data objects are not written out to the final `yoda` file. This option writes out such objects.

## Do it yourself!

If you need features which are not implemented here, or have your own specific mergign requirements, 
then you could also make a local copy of `yodamerge` and use the very powerful `Python` API to implement
your own merging prescriptions. 

