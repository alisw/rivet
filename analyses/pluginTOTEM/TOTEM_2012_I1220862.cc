// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// TOTEM elastic and total cross-section measurement
  class TOTEM_2012_I1220862 : public Analysis {
  public:

    TOTEM_2012_I1220862()
      : Analysis("TOTEM_2012_I1220862")
    {    }


    void init() {
      declare(ChargedFinalState(), "CFS");
      book(_hist_tlow  ,1, 1, 1);
      book(_hist_thigh ,2, 1, 1);
      book(_hist_sigma ,3, 1, 1);
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() > 2) MSG_DEBUG("Final state includes more than two charged particles!");
      _hist_sigma->fill(sqrtS()/GeV);

      for (const Particle& p : cfs.particles(Cuts::eta > 0)) { // && Cuts::pid == PID::PROTON)) {
        if (p.pid() != PID::PROTON) continue;
        const double t = sqr(p.pT());
        _hist_tlow->fill(t);
        _hist_thigh->fill(t);
      }
    }


    void finalize() {
      normalize(_hist_tlow, crossSection()/millibarn);
      normalize(_hist_thigh, crossSection()/millibarn);
      normalize(_hist_sigma, crossSection()/millibarn);
    }


  private:

    Histo1DPtr _hist_tlow, _hist_thigh, _hist_sigma;

  };


  DECLARE_ALIASED_RIVET_PLUGIN(TOTEM_2012_I1220862, TOTEM_2012_002);

}
