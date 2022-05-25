// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TOTEM_2012_I1115294 : public Analysis {
  public:

    TOTEM_2012_I1115294()
      : Analysis("TOTEM_2012_I1115294")
    {    }


  public:

    void init() {
      ChargedFinalState cfsm((Cuts::etaIn(-6.50, -5.35) && Cuts::pT >=  40.*MeV));
      ChargedFinalState cfsp((Cuts::etaIn( 5.35,  6.50) && Cuts::pT >=  40.*MeV));
      declare(cfsm, "CFSM");
      declare(cfsp, "CFSP");

      book(_h_eta ,1, 1, 1);
      book(_sumofweights, "sumofweights");
    }


    void analyze(const Event& event) {
      const ChargedFinalState cfsm = apply<ChargedFinalState>(event, "CFSM");
      const ChargedFinalState cfsp = apply<ChargedFinalState>(event, "CFSP");

      if (cfsm.size() == 0 && cfsp.size() == 0) vetoEvent;

      _sumofweights->fill();

      for (const Particle& p : cfsm.particles() + cfsp.particles()) {
        _h_eta->fill(p.abseta());
      }

    }


    void finalize() {
      scale(_h_eta, 0.5 / *_sumofweights);
    }


  private:

    CounterPtr _sumofweights;
    Histo1DPtr _h_eta;


  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TOTEM_2012_I1115294);

}
