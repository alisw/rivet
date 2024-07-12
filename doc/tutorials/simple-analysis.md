## Instructions for writing a Rivet routine/analysis

At this point you have probably run Rivet a few times, maybe different analyses, different generators, and played around with plotting options. To really utilize the full power of the framework, one does, however, need to write analysis code -- it is an analysis framework after all.

The best tip for writing analyses, is to find an existing similar analysis, from the large library of already existing ones, and take inspiration from that. But even then, it is important to have the basics right. 

### Writing the analysis code

Here we are going to write a new analysis for use with Rivet. This is done "stand-alone", i.e. you don't have to modify the code of Rivet itself: in fact, you can follow these instructions using a system install of Rivet to which you have no write permissions.

All analysis routines are implemented as sub-classes of the Rivet "Analysis" class: pretty much all the magic that binds the analysis object into the Rivet system is handled in this base class, meaning that your code can really concentrate on implementing the physics goals of the analysis. 


### The analysis "wizard"

You could make your analysis by copying some example code and then going through a load of search and replace fiddling, but in fact there is a much easier way: the `rivet-mkanalysis` script. This script is installed along with the rest of the Rivet system, and will generate template analysis files for you.

You can get some help info by running `rivet-mkanalysis --help`, but the basic usage (to generate the files in your current directory) is `rivet-mkanalysis MY_ANALYSIS_NAME`. A three part name, separated by underscores, is a Rivet convention that we recommend you to use: the first part is the experiment name, the second is the year of publication, and the third is the ID code for the corresponding paper in the [Inspire HEP database](http://inspirehep.net), preceded by an "I". You can get the Inspire ID from an Inspire record page by looking at the URL: it will be the numerical trailing part of the address following `record/`, e.g. 849050 in the record [http://inspirehep.net/record/849050.](http://inspirehep.net/record/849050)

So, for example, `ATLAS_2010_I849050` would be the standard name for the analysis link given above... although in fact that one (the first ATLAS minimum bias paper) is in the Rivet collection under the name `ATLAS_2010_S8591806`, which uses the older ''SPIRES'' database ID. Please use Inspire IDs rather than SPIRES ones for new analyses -- we intend to update all analyses to the Inspire naming in a future release.

Running the `rivet-mkanalysis` script with the appropriate analysis name will have generated a `.cc` C++ source file template, and template metadata files for information about the analysis (`.info`) and specifications of titles, axis labels, etc. for the plots which the analysis will produce (`.plot`). These templates will include, if possible, extra analysis metadata such as a BibTeX publication entry in the `.info` file.

### Structure

For simplicity, Rivet analysis classes are usually written in just one `.cc` file, i.e. no header declaration. This is because classes are almost always not inherited from, and all that the Rivet system needs to know is that it can be treated as an `Analysis*` pointer: avoiding header files makes everything more compact and removes a source of errors and annoyance.

An analysis has the following components:
 * a no-argument constructor;
 * three analysis event loop methods: `init`, `analyze` and `finalize`;
 * a minimal hook into the plugin system

It is also possible to add some metadata methods which describe the analysis, references to publications, experiment, etc., but we strongly recommend that you put this information into the "YAML" format (see http://www.yaml.org) `.info` template that the `rivet-mkanalysis` script generated for you instead: this way the code will remain clean and minimal, and you can update the metadata without needing to recompile. All analyses bundled with Rivet store their metadata in external files.

Useful analyses also contain member variables for the analysis: event weight counters and histograms are the most common of these. Conventionally, we declare the class member variables with a leading underscore: see the [Coding Style Guide](codingstyle.md) for more information on our recommended uniform coding style. Histogram pointer members (for which we use special smart pointers with clever machinery inside) are preferred to start with an "h", e.g. `Histo1DPtr _h_pT`.

### Implementation

The constructor and three event loop methods are used for the following:

 * Constructor: set whether the generator cross-section is needed. Minimal!
 * `init`: book histograms, initialise counters, etc.
 * `analyze`: select particles, filter according to cuts, loop over combinations, construct observables, fill histograms. This is where the per-event aspect of the analysis algorithm goes.
 * `finalize`: normalize/scale/divide histograms, tuples, etc.

This probably looks similar to every analysis system you've ever used, so hopefully you're not worried about Rivet being weird or difficult to learn ;-)

Rivet provides implementations of many calculational tools, called "projections". These are just observable calculator objects with a silly name, so don't get worried. (They automatically cache their results, to make Rivet automatically efficient, but you don't have to worry about that since it's, well, automatic.) The projections are used by calling the analysis' `apply(event)` method. This will return a const reference to the completed projection object and takes the type of the reference as a template argument, e.g.
```
  const FinalState& cfs = apply<FinalState>(event, "Tracks");
```
The name "Tracks" here will have been registered in the `init` method as referring to a projection of type "ChargedFinalState" --- a calculator which provides a list of charged particles with certain basic cuts applied. This is done via the `declare` method. Note that a) you don't have to manage the memory yourself, and b) polymorphism via the reference is both allowed and encouraged. If b) means nothing to you, don't worry... we just want to reassure C++ fiends who might think we're cramping their style!

### Example

Here is an example of the whole Rivet analysis shebang. We've compressed it into a single .cc file since the `analyze` method is nice and short and there is no reason to make a header:

```
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  
  class MyAnalysis : public Analysis {
  public:
    
    /// Default constructor
    MyAnalysis() : Analysis("MYANALYSIS") {   }
    
    
    /// @name Analysis methods
    //@{
    void init() {
      const FinalState fs(Cuts::abseta < 5);
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");
      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV), "Tracks");
    }
    
    void analyze(const Event& event) {
      const Jets& jets = apply<ChargedFinalState>(event, "Jets")
        .jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 4.4);
      MSG_DEBUG("Jet multiplicity = " << jets.size());

      const Particles& trks = apply<FinalState>(event, "Tracks").particles();
      MSG_DEBUG("Track multiplicity = " << trks.size());
    }

    // No histos, so no need for a finalize()!

    //@}
    
  };

  // Magic required by the plugin system 
  DECLARE_RIVET_PLUGIN(MyAnalysis);
  
}
```

### Cut objects

Note the use of objects in the `Cuts` namespace to specify kinematic cuts on particles or jets selected by projections, or returned from them as lists. These predefined objects of type `Rivet::Cut` can be combined together using arbitrary combinations of logical operators, with the combined object also being of type `Cut`.

Many functions in Rivet accept a (potentially compound) `Cut` as an argument, so this is a very flexible, unambiguous, and human-readable way to express analysis selection logic. There is not much to know from the user point of view beyond what you see above!

The standard Rivet predefined cuts are (all in the `Rivet::Cuts` namespace): `pT`, `Et`, `mass`, `phi`, `eta`, `abseta`, `rap`, `absrap`.


### Compiling and linking
 
To use your new analysis, you need to build it into a Rivet analysis plugin library, with a name of the form `Rivet*.so` library. You can do this manually, but to make life easier there is again a helper script, used as follows:
`rivet-buildplugin RivetMyAnalyses.so MyAnalysis.cc MyOtherAnalysis.cc # etc.`

Note that the name of the library has to start with the word "Rivet" or it will not get loaded at runtime. By default, if no ".so" first argument is given, the name =RivetAnalysis.so= will be used.


### Running

You can now use your new analysis right away. Provided that the `RivetMyAnalysis.so` shared library file, or a similarly-named symbolic link to it, is in a directory listed in your `RIVET_ANALYSIS_PATH` environment variable, it will work right away with the `rivet` command:
```
  > ls
  RivetMyAnalysis.so MyAnalysis.cc
  > export RIVET_ANALYSIS_PATH=$PWD
  > rivet --list-analyses
  [...]
  MYANALYSIS
  > rivet --show-analysis MYANALYSIS
  MyAnalysis
  ==========

  Spires ID: NONE
  Spires URL: http://www.slac.stanford.edu/spires/find/hep/www?rawcmd=key+NONE
  Experiment: NONE
  Year of publication: NONE

  Description:
    A do-nothing analysis for demonstrating how to make a plugin

  References:
```

Alternatively, you can use the `--analysis-path` flag to `rivet`:
```
  > rivet --list-analyses --analysis-path=$PWD
```

### Making it useful

Hopefully that's enough to get you started. The other main things to learn are booking (and "auto-booking") of histograms and other data objects, and use of the Rivet projections and analysis objects. For this, we recommend that you take a look at the code of some of the standard analyses, and read more information about projections and histogramming on this wiki and in the Rivet PDF manual.
