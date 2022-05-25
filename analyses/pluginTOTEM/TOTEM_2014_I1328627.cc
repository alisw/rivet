// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TOTEM_2014_I1328627 : public Analysis {
  public:


    TOTEM_2014_I1328627()
      : Analysis("TOTEM_2014_I1328627")
    {    }



    void init() {
      ChargedFinalState cfsm(Cuts::etaIn(-7.0, -6.0));
      ChargedFinalState cfsp(Cuts::etaIn( 3.7,  4.8));
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
      scale(_h_eta, 1./ *_sumofweights);
    }


  private:

    CounterPtr _sumofweights;
    Histo1DPtr _h_eta;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TOTEM_2014_I1328627);


}
