// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"


namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_2021_I1972986 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1972986);


    /// Book histograms and initialize projections:
    void init() {

      // Initialize the projections
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");

      // Book sets of histograms, binned in absolute rapidity
      Histo1DPtr tmp;
      for(int y = 0; y < 4; ++y) {
         _hist_sigmaAK4.add(0.5*y, 0.5*(y+1), book(tmp,y+1, 1, 1)); // d0?-x01-y01
         _hist_sigmaAK7.add(0.5*y, 0.5*(y+1), book(tmp,20+y+1, 1, 1)); // d2?-x01-y01
      }
    }


    /// Per-event analysis
    void analyze(const Event &event) {

      // AK4 jets
      const FastJets& fjAK4 = apply<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jets(Cuts::ptIn(97*GeV, 3103*GeV) && Cuts::absrap < 2.0);
      for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4.fill(j.absrap(), j.pT());
      }

      // AK7 jets
      const FastJets& fjAK7 = apply<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jets(Cuts::ptIn(97*GeV, 3103*GeV) && Cuts::absrap < 2.0);
      for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7.fill(j.absrap(), j.pT());
      }

    }


    // Finalize
    void finalize() {
      /// @note Extra division factor is the *signed* dy, i.e. 2*d|y|
      _hist_sigmaAK4.scale(crossSection()/picobarn/sumOfWeights()/2.0, this);
      _hist_sigmaAK7.scale(crossSection()/picobarn/sumOfWeights()/2.0, this);
    }


    /// @name Histograms
    //@{
    BinnedHistogram _hist_sigmaAK4;
    BinnedHistogram _hist_sigmaAK7;
    //@}

  };


  // This global object acts as a hook for the plugin system.
  RIVET_DECLARE_PLUGIN(CMS_2021_I1972986);

}
