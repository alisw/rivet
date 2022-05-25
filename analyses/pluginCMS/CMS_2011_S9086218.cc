// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  // Inclusive jet pT
  class CMS_2011_S9086218 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S9086218);


    /// @{

    /// Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      // Initialize the projectors:
      declare(FastJets(fs, FastJets::ANTIKT, 0.5),"Jets");

      // Book histograms:
      {Histo1DPtr tmp; _hist_sigma.add(0.0, 0.5, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(0.5, 1.0, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(1.0, 1.5, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(1.5, 2.0, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(2.0, 2.5, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(2.5, 3.0, book(tmp, 6, 1, 1));}
    }

    /// Analysis
    void analyze(const Event &event) {
      const double weight = 1.0;
      const FastJets& fj = apply<FastJets>(event,"Jets");
      const Jets& jets = fj.jets(Cuts::ptIn(18*GeV, 1100.0*GeV) && Cuts::absrap < 4.7);

      // Fill the relevant histograms:
      for(const Jet& j : jets) {
        _hist_sigma.fill(j.absrap(), j.pT(), weight);
      }
    }

    /// Finalize
    void finalize() {
      _hist_sigma.scale(crossSection()/sumOfWeights()/2.0, this);
    }

    /// @}


  private:

    BinnedHistogram _hist_sigma;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S9086218, CMS_2011_I902309);

}
