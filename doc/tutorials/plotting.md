## Plotting

Rivet comes with three commands --- `rivet-mkhtml`, `compare-histos` and `make-plots` --- for comparing and plotting data files. These commands produce nice comparison plots of publication quality from the YODA format text files.

The high-level program `rivet-mkhtml` will automatically create a plot
webpage from the given YODA files. It searches for reference data automatically
and uses the other two commands internally. Example:

```
rivet-mkhtml withUE.yoda:'Title=With UE' withoutUE.yoda:'LineColor=blue'
```

Run `rivet-mkhtml --help` to find out about all features and options.

Plotting options will be taken from `*.plot` files installed in Rivet's share directory.
These contain plotting instructions according to the `make-plots` syntax as documented in
[http://rivet.hepforge.org/make-plots.html.](http://rivet.hepforge.org/make-plots.html)
Such files can also be written for any plugin analysis and will be found if they are in the `RIVET_ANALYSIS_PATH`.

You can also run the other two commands separately:

`compare-histos` will accept a number of YODA files as input (ending in
`.yoda`, identify which plots are available in them, and combine the MC
and reference plots appropriately into a set of plot data files ending with
`.dat`. More options are described by running `compare-histos --help`.

Incidentally, the reference files for each Rivet analysis are to be found
in the installed Rivet shared data directory, `installdir/share/Rivet`.
You can find the location of this by using the `rivet-config` command:

```
rivet-config --datadir
```

You can plot the created data files using the `make-plots` command:

```
make-plots --pdf *.dat
```

The `--pdf` flag makes the output plots in PDF format: by default the output
is in postscript (`.ps`), and flags for conversion to EPS and PNG are also
available.

Note that the plotting tools internally use LaTeX for drawing, and for very
complex plots it might sometimes fail with an error message like
"TeX memory exceeded" (or "DVI file can't be opened"). In such
a case it is recommended to increase the allowed TeX memory size as described
e.g. in the [pgfplots manual](http://pgfplots.sourceforge.net/pgfplots.pdf)
in Section 6.2.

Furthermore, imagemagick's `convert` is used to convert .pdf files into .png
files for the plot website. Some versions of imagemagick might cause error
messages like `not authorized @ error/constitute.c/ReadImage/412`, which can
[be fixed](https://stackoverflow.com/q/52998331/3094872) by giving `rights="read|write"` for both PS and PDF formats in
`/etc/ImageMagick-7/policy.xml`.
