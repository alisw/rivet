# Running on a subset of available multiweights 

By default Rivet will run over all the available multiweights in can find in the input file.
It might be beneficial to only run on a subset of the available multiweights, especially
when the output file would be rather large and the weights aren't actually needed.
The following flags can help:

* passing the `--skip-multiweights` (or its alias`--skip-weights`) flag to Rivet will skip
all multiweights apart from the nominal.
* the `--match-weights` flag can be used to select a subset of the weights, 
e.g. `--match-weights=WeightName1,WeightName2`, or using a suitable regular expressions,
e.g. `--match-weights=MUR.*MUF.*PDF123456`. Note that the default weight can never be de-selected.
* the `--unmatch-weights` flag can be used to de-select a subset of the weights, 
e.g. `--match-weights=WeightName1,WeightName2`, or using a suitable regular expressions,
e.g. `--match-weights=MUR1_MUF1_PDF123.*`. Note that the default weight can never be de-selected.

When both flags are passed, `--match-weights` is evaluated first, and `--unmatch-weights` is applied
to the "surviving weights". Whenever `--skip-multiweights` is passed, it takes precedence 
(i.e. `--match-weights` and `--unmatch-weight` would be ignored in that case).

