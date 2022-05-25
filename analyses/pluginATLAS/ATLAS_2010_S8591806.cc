// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief ATLAS minimum bias analysis at 900 GeV
  /// @author Frank Siegert
  class ATLAS_2010_S8591806 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2010_S8591806);


    void init() {
      ChargedFinalState cfs((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  0.5*GeV));
      declare(cfs, "CFS");

      book(_h_dNch_deta ,2, 1, 1);
      book(_h_dNch_dpT ,3, 1, 1);
      book(_h_dNevt_dNch ,4, 1, 1);
      book(_p_meanpT_Nch ,5, 1, 1);

      book(_Nevt_after_cuts, "nevt_pass");
    }


    void analyze(const Event& event) {
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      if (charged.size() < 1) {
        vetoEvent;
      }
      _Nevt_after_cuts->fill();

      _h_dNevt_dNch->fill(charged.size());
      for (const Particle& p : charged.particles()) {
        double pT = p.pT()/GeV;
        _h_dNch_deta->fill(p.eta());
        _h_dNch_dpT->fill(pT, 1.0/pT);
        _p_meanpT_Nch->fill(charged.size(), pT);
      }
    }


    void finalize() {
      double deta = 5.0;
      scale(_h_dNch_deta, 1.0/_Nevt_after_cuts->val());
      scale(_h_dNch_dpT, 1.0/_Nevt_after_cuts->val()/TWOPI/deta);
      scale(_h_dNevt_dNch, 1.0/_Nevt_after_cuts->val());
    }


  private:

    Histo1DPtr _h_dNch_deta;
    Histo1DPtr _h_dNch_dpT;
    Histo1DPtr _h_dNevt_dNch;
    Profile1DPtr  _p_meanpT_Nch;

    CounterPtr _Nevt_after_cuts;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2010_S8591806, ATLAS_2010_I849050);

}
