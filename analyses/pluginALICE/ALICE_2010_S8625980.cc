// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2010_S8625980 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2010_S8625980);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs((Cuts::etaIn(-1.0, 1.0)));
      declare(cfs, "CFS");

      if (isCompatibleWithSqrtS(900)) {
        book(_h_dN_deta    ,4, 1, 1);
      } else if (isCompatibleWithSqrtS(2360)) {
        book(_h_dN_deta    ,5, 1, 1);
      } else if (isCompatibleWithSqrtS(7000)) {
        book(_h_dN_deta    ,6, 1, 1);
        book(_h_dN_dNch    ,3, 1, 1);
      }
      book(_Nevt_after_cuts, "Nevt_after_cuts");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      if (charged.size() < 1) {
        vetoEvent;
      }
      _Nevt_after_cuts->fill();


      for (const Particle& p : charged.particles()) {
        const double eta = p.eta();
        _h_dN_deta->fill(eta);
      }

      if (isCompatibleWithSqrtS(7000)) {
        _h_dN_dNch->fill(charged.size());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      if (isCompatibleWithSqrtS(7000)) {
        normalize(_h_dN_dNch);
      }
      scale(_h_dN_deta, 1.0/ *_Nevt_after_cuts);

    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_dN_deta;
    Histo1DPtr _h_dN_dNch;
    CounterPtr _Nevt_after_cuts;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALICE_2010_S8625980, ALICE_2010_I852264);

}
