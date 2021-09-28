// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  class CMSTOTEM_2014_I1294140 : public Analysis {
  public:

    CMSTOTEM_2014_I1294140()
      : Analysis("CMSTOTEM_2014_I1294140")
    {     }


    void init() {
      ChargedFinalState cfs((Cuts::etaIn(-7.0, 7.0)));
      declare(cfs, "CFS");

      book(_Nevt_after_cuts_or, "Nevt_or");
      book(_Nevt_after_cuts_and, "Nevt_and");
      book(_Nevt_after_cuts_xor, "Nevt_xor");

      if (fuzzyEquals(sqrtS(), 8000*GeV, 1E-3)) {
        book(_h_dNch_dEta_OR ,1, 1, 1);
        book(_h_dNch_dEta_AND ,2, 1, 1);
        book(_h_dNch_dEta_XOR ,3, 1, 1);
      }
    }


    void analyze(const Event& event) {
      // Count forward and backward charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      int count_plus = 0, count_minus = 0;
      for (const Particle& p : charged.particles()) {
        if (inRange(p.eta(),  5.3,  6.5)) count_plus++;
        if (inRange(p.eta(), -6.5, -5.3)) count_minus++;
      }

      // Cut combinations
      const bool cutsor  = (count_plus > 0 || count_minus > 0);
      const bool cutsand = (count_plus > 0 && count_minus > 0);
      const bool cutsxor = ( (count_plus > 0 && count_minus == 0) || (count_plus == 0 && count_minus > 0) );

      // Increment counters and fill histos
      if (cutsor)  _Nevt_after_cuts_or  ->fill();
      if (cutsand) _Nevt_after_cuts_and ->fill();
      if (cutsxor) _Nevt_after_cuts_xor ->fill();
      for (const Particle& p : charged.particles()) {
        if (cutsor)  _h_dNch_dEta_OR ->fill(p.abseta());
        if (cutsand) _h_dNch_dEta_AND->fill(p.abseta());
        if (cutsxor) _h_dNch_dEta_XOR->fill(p.abseta());
      }

    }


    void finalize() {
      scale(_h_dNch_dEta_OR,  0.5 / *_Nevt_after_cuts_or);
      scale(_h_dNch_dEta_AND, 0.5 / *_Nevt_after_cuts_and);
      scale(_h_dNch_dEta_XOR, 0.5 / *_Nevt_after_cuts_xor);
    }


  private:

    Histo1DPtr _h_dNch_dEta_OR, _h_dNch_dEta_AND, _h_dNch_dEta_XOR;
    CounterPtr _Nevt_after_cuts_or, _Nevt_after_cuts_and, _Nevt_after_cuts_xor;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMSTOTEM_2014_I1294140);

}
