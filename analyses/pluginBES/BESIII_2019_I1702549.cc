// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief q^2 in D_s+ -> K0 e+ nu_e
  class BESIII_2019_I1702549 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2019_I1702549);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2, 1, 1, 1);
      book(_nD,"/TMP/nD");
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over D+ mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==431 )) {
        _nD->fill();
        if(p.pid()==431 && isSemileptonicDecay(p, {PID::K0, PID::EPLUS, PID::NU_E}) ) {
      	  _h_q2->fill(q2(p, PID::K0));
        }
        else if(p.pid()==431 && isSemileptonicDecay(p, {130, PID::EPLUS, PID::NU_E}) ) {
      	  _h_q2->fill(q2(p, 130));
        }
        else if(p.pid()==431 && isSemileptonicDecay(p, {310, PID::EPLUS, PID::NU_E}) ) {
      	  _h_q2->fill(q2(p, 310));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       // normalise to width in inverse ns
       scale(_h_q2, 1./0.504e-3/ *_nD);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_q2;
    CounterPtr _nD;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BESIII_2019_I1702549);

}
