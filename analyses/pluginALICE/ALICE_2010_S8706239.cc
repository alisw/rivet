// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2010_S8706239 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ALICE_2010_S8706239()
      : Analysis("ALICE_2010_S8706239")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs((Cuts::etaIn(-0.8, 0.8) && Cuts::pT >=  0.15));
      declare(cfs, "CFS");

      book(_h_pT ,4, 1, 1);

      book(_h_pT_Nch_015 ,11, 1, 1);
      book(_h_pT_Nch_05  ,12, 1, 1);

      book(_Nevt_after_cuts,"Nevt_after_cuts");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      _Nevt_after_cuts->fill();

      // Get number of particles that fulfill certain pT requirements
      int Nch_015 = 0;
      int Nch_05  = 0;
      for (const Particle& p : charged.particles()) {
        double pT = p.pT()/GeV;
        if (pT < 4.0) Nch_015++;
        if (pT > 0.5  && pT < 4.0) Nch_05++;
      }

      // Now we can fill histograms
      for (const Particle& p : charged.particles()) {
        double pT = p.pT()/GeV;
        if (pT < 4.0) _h_pT_Nch_015 ->fill(Nch_015, pT);
        if (pT > 0.5  && pT < 4.0) _h_pT_Nch_05  ->fill(Nch_05,  pT);

      // To get the Yield, fill appropriate weight 1/(2PI * pT * d eta)
        _h_pT->fill(pT, 1.0 /(TWOPI*pT*1.6) );
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pT, 1.0/ *_Nevt_after_cuts);
    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_pT;

    Profile1DPtr _h_pT_Nch_015 ;
    Profile1DPtr _h_pT_Nch_05  ;

    CounterPtr _Nevt_after_cuts;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2010_S8706239);

}
