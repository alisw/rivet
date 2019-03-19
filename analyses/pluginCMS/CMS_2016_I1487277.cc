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
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.7),"Jets");

      // Book histograms:

      
      _hist_sigma.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _hist_sigma.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 1));
      _hist_sigma.addHistogram(1.0, 1.5, bookHisto1D(3, 1, 1));
      _hist_sigma.addHistogram(1.5, 2.0, bookHisto1D(4, 1, 1));
      _hist_sigma.addHistogram(2.0, 2.5, bookHisto1D(5, 1, 1));
      _hist_sigma.addHistogram(2.5, 3.0, bookHisto1D(6, 1, 1));
      _hist_sigma.addHistogram(3.2, 4.7, bookHisto1D(7, 1, 1));

    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = event.weight();
      const FastJets &fj = applyProjection<FastJets>(event,"Jets");
      const Jets& jets = fj.jets(Cuts::ptIn(18*GeV, 5000.0*GeV) && Cuts::absrap < 5.2);

      // Fill the relevant histograms:
      foreach(const Jet &j, jets) {
        _hist_sigma.fill(j.absrap(), j.pT(), weight);
      }
    }

    // Finalize
    void finalize() {
    _hist_sigma.scale(crossSection()/sumOfWeights()/2.0, this);
    }

  private:
    BinnedHistogram<double> _hist_sigma;
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
