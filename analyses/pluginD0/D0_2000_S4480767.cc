// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  class D0_2000_S4480767 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2000_S4480767);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      WFinder wf(fs, Cuts::abseta < 5, PID::ELECTRON, 0.0*GeV, 200.0*GeV, 0.0*GeV, 0.2);
      declare(wf, "WFinder");

      book(_h_W_pT ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const WFinder& wf = apply<WFinder>(event, "WFinder");
      if (wf.bosons().size() == 0) vetoEvent;

      _h_W_pT->fill(wf.bosons()[0].pT()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_W_pT, crossSection()/sumOfWeights());
    }

    /// @}


  private:

    /// Histogram
    Histo1DPtr _h_W_pT;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2000_S4480767, D0_2000_I535017);

}
