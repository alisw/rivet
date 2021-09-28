// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // Inclusive jet pT
  class CMS_2016_I1487277 : public Analysis {
  public:

    // Constructor
    CMS_2016_I1487277() : Analysis("CMS_2016_I1487277") {}


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      // Initialize the projectors:
      declare(FastJets(fs, FastJets::ANTIKT, 0.7),"Jets");

      // Book histograms:

      
      {Histo1DPtr tmp; _hist_sigma.add(0.0, 0.5, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(0.5, 1.0, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(1.0, 1.5, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(1.5, 2.0, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(2.0, 2.5, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(2.5, 3.0, book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _hist_sigma.add(3.2, 4.7, book(tmp, 7, 1, 1));}

    }

    // Analysis
    void analyze(const Event &event) {
      const FastJets &fj = applyProjection<FastJets>(event,"Jets");
      const Jets& jets = fj.jets(Cuts::ptIn(18*GeV, 5000.0*GeV) && Cuts::absrap < 5.2);

      // Fill the relevant histograms:
      for(const Jet &j : jets) {
        _hist_sigma.fill(j.absrap(), j.pT());
      }
    }

    // Finalize
    void finalize() {
    _hist_sigma.scale(crossSection()/sumOfWeights()/2.0, this);
    }

  private:
    BinnedHistogram _hist_sigma;
    Histo1DPtr _hist_ptbins_y1;
    Histo1DPtr _hist_ptbins_y2;
    Histo1DPtr _hist_ptbins_y3;
    Histo1DPtr _hist_ptbins_y4;
    Histo1DPtr _hist_ptbins_y5;
    Histo1DPtr _hist_ptbins_y6;
    Histo1DPtr _hist_ptbins_y7;

  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2016_I1487277);

}
