## Histogramming in Rivet

Histogramming in Rivet is currently handled via the YODA data analysis interfaces: see (https://yoda.hepforge.org/)[https://yoda.hepforge.org] for more information on YODA.

### Booking histograms

Most of the time you will want to book histograms from within an Analysis. Rivet provides machinery to handle the both contruction of the histogram and making the histogram available for writing. There are several ways of booking a histogram:
 
 * book evenly spaced bins in a range, by specifying the endpoints and the number of bins
 * pass a specific `std::vector` of bin edges
 * automatically book with the right bin edges based on the reference data file.

Any "official" Rivet analyses will use the latter method, and we recommend that you do so, too. If you are starting a new analysis, and the reference data is found in HepData (the INSPIRE record will provide a link to "Data: HepData"), it can be downloaded from there -- and will even be done automatically if an analysis skeleton is written by the `rivet-mkanalysis` script.

The `book` method takes a pointer to any histogram type (all derived from `YODA::AnalysisObject) as the first argument, and will initialize the object behind the pointer accordingly.

### Tell me more about this auto-booking thing...

Maybe this isn't so obvious after all! The idea is that most Rivet analyses should be comparable with experimental data, such as that in the HepData database. 

Since the MC and ref data must have the same binnings to be meaningfully compared, and since encoding long lists of bin edges into your code is annoying, error-prone and ugly, our booking system will use the reference data files as a template from which to book the MC histogram. The reference files will be searched for in the installation path of Rivet. The internal path to the reference histograms is the same as for the MC histograms, but all inside a top-level virtual directory called "REF". For example, MC histo `/MY_ANALYSIS/my-histo` can be auto-booked from reference histo `/REF/MY_ANALYSIS/my-hist`

For those histograms dumped from HepData, the histogram naming system is a little strange: Every 1D histogram that you might want to make is stored in HepData as a combination of an x-axis and a y-axis for a given dataset --- HepData datasets can have multiple axes of either kind, so that e.g. several measurements binned the same way can be encoded as multiple y-axes on a single x-axis. So, armed with the dataset, x-axis and y-axis IDs, if you call:
```
// In MyAnalysis class definition:
Histo1DPtr myHisto;

// In implementation:
void MyAnalysis::init() {
  ...
  book(myHisto, 1, 1, 2);
  ...
}
```
then you'll get `myHisto` initialized with the right binning. 

### Filling histograms

Histograms are usually filled in your analysis' `analyze()` method, using the normal YODA fill() method, e.g.
```
histMyDistribution->fill(value);
```
Note that contrary to what you may be used to from ROOT or other analysis frameworks, it is not neccessary to fill
the event weight into the histogram. This is handled automatically by Rivet.

