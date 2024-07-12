// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2010_S8624100 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2010_S8624100);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs05((Cuts::etaIn(-0.5, 0.5)));
      ChargedFinalState cfs10((Cuts::etaIn(-1.0, 1.0)));
      ChargedFinalState cfs13((Cuts::etaIn(-1.3, 1.3)));
      declare(cfs05, "CFS05");
      declare(cfs10, "CFS10");
      declare(cfs13, "CFS13");

      if (isCompatibleWithSqrtS(900)) {
        book(_h_dN_dNch_05    ,11, 1, 1);
        book(_h_dN_dNch_10    ,12, 1, 1);
        book(_h_dN_dNch_13    ,13, 1, 1);
      } else if (isCompatibleWithSqrtS(2360)) {
        book(_h_dN_dNch_05    ,17, 1, 1);
        book(_h_dN_dNch_10    ,18, 1, 1);
        book(_h_dN_dNch_13    ,19, 1, 1);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
       const ChargedFinalState& charged_05 = apply<ChargedFinalState>(event, "CFS05");
      const ChargedFinalState& charged_10 = apply<ChargedFinalState>(event, "CFS10");
      const ChargedFinalState& charged_13 = apply<ChargedFinalState>(event, "CFS13");

      _h_dN_dNch_05->fill(charged_05.size());
      _h_dN_dNch_10->fill(charged_10.size());
      _h_dN_dNch_13->fill(charged_13.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_dN_dNch_05);
      normalize(_h_dN_dNch_10);
      normalize(_h_dN_dNch_13);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_dN_dNch_05;
    Histo1DPtr _h_dN_dNch_10;
    Histo1DPtr _h_dN_dNch_13;
    /// @}


  };


  RIVET_DECLARE_ALIASED_PLUGIN(ALICE_2010_S8624100, ALICE_2010_I852450);

}
