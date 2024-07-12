# Fun with multiweights when plotting YODA files (work in progress)

## Weight names

The weight name should be a reasonably concise string - please don't write an essay. 
It should be a sequence of alphanumeric characters and can include the following special characters:
`=`, `_`, `.`, `+` and `-`. Other special characters are currently not supported.

White space characters are in principle supported but have been a source of endless headaches 
for Rivet users and Monte Carlo users in general: please don't use white space in weight names.

Ideally, the weight name should include the LHAPDF ID where applicable, e.g. `MUR2.0_MUF1.0_PDF123456`.

## Plotting a subset of the weights

In order to select a subset of the weight names for plotting, specifiy the desired subset 
using a comma-separated list as follows:
```
rivet-mkhtml --errs file.yoda:"Variations=Weight1,Weight2,MUR.*MUF.*PDF123456"
```
Note that regular expressions are also supported.
Unless there is a multiweight called "none" in the file, 
`Variations=none` could be passed to veto all multiweights
for a specific file. Passing the `--no-weights` flag will
remove the multiweights for all input files.


## Combining the weights into uncertainty bands

In general, combining weights will always require input from the user at some level, 
especially when a more complicated prescription for how to combine the weights is required. 
Some basic/frequently used options are available in the scripts, but as always — more complicated 
things will require the user to write their own scripts e.g. using the excellent YODA Python API 
to combine the histograms prior to the plotting step.

For simple things, two options are available which can be used several times for a given curve.
All band uncertainty components will currently be added in quadrature together with the 
MC statistical uncertainty in order to produce a single "total uncertainty band" for the given curve.
The options are:

### Envelopes

The `BandComponentEnv` will take all the variations whose weight name matches the 
comma-separated list (or regular expression) and work out the envelope from them.

```
rivet-mkhtml --errs file1.yoda:"BandComponentEnv=Weight1,Weight2" file2.yoda:"BandComponentEnv=MUR.*MUF.*PDF123456"
```

If multiple envelopes are required, just leave a space between the groups of weight names,
for which separate envelopes will be worked out (and then summed in quadrature).:

```
rivet-mkhtml --errs file.yoda:"BandComponentEnv=Weight1,Weight2 MUR.*MUF.*PDF123456"
```

### PDF bands

In analogy to the envelopes, there is also a `BandComponentPDF` tag that follows a similar structure.

```
rivet-mkhtml --errs file.yoda:"BandComponentPDF=MUR1_MUF1_PDF123.*"
```
However, this will only work if 
* `LHAPDF` is installed such that it can import the `LHAPDF` Python API, 
* the PDF used is known to the `LHAPDF` tool, 
* the exact number of weights names associated with the PDF set are captured and 
* the `LHAPDF` IDs are actually specified in the weight name starting with "PDF" like in the example above.


## Choosing a (different) central value

Rivet tries hard to identify which one of the weights is the default weight representing the central value.
However, in the absence of a HEP-wide standard for weight names, this might sometimes fail.
In that case, Rivet itself will happily process the input file and just treat all weights as variation weights
(i.e. it will also attach the usual `[WeightName]` to the output histogram name) and it's up to the user
to let the plotting script know which weight should be treated as the central value.
This is achieved as follows:

```
rivet-mkhtml --errs file.yoda:"DefaultWeight=WeightName"
```

Note that the same trick can be used to choose a different central value at plotting time
even if Rivet managed to identify the actual default weight at run time.

## Comparing curves within the same file

Sometimes you may just want to compare the nominal to a variation weight from the same file.
It is possible to trick the plotting tools into thinking they are looking at different files
by passing the same file twice but giving it a different `Name`, e.g. like so

```
rivet-mkhtml --errs file.yoda:"Name=nominal:Variations=none" file.yoda:"Name=variation:DefaultWeight=Weight1:Variations=none"
```

