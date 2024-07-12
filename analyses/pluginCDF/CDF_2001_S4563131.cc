// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run I inclusive jet cross-section
  class CDF_2001_S4563131 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2001_S4563131);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");
      book(_h_ET ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 40*GeV && Cuts::abseta >= 0.1 && Cuts::abseta <= 0.7, cmpMomByEt);
      for (const Jet& jet : jets) {
        _h_ET->fill(jet.Et());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double deta = 1.2;
      scale(_h_ET, crossSection()/sumOfWeights()/deta/nanobarn);
    }

    //@}


  private:

    /// Histogram
    Histo1DPtr _h_ET;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2001_S4563131, CDF_2001_I552797);

}
