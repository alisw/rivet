// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class MY_ANALYSIS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MY_ANALYSIS);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      _lmode = 0; // default accepts either channel
      if ( getOption("LMODE") == "EL" )  _lmode = 1;
      if ( getOption("LMODE") == "MU" )  _lmode = 2;

      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      book(_h["mll"], "mass_ll", mll_bins);
      //book(_h["jets_excl"],  "jets_excl",   6, -0.5,  5.5);
      //book(_h["bjets_excl"], "bjets_excl",  3, -0.5,  2.5);
      //book(_h["HT"],         "HT",          6,  20., 110.);
      //book(_h["pTmiss"],     "pTmiss",     10,   0., 100.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// Todo: Reconstruct the dilepton invariant mass to fill the histogram
      // ... 
      _h["mll"]->fill(1.0);

    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection() / sumOfWeights();
      scale(_h, sf);

    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    size_t _lmode;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MY_ANALYSIS);
}
