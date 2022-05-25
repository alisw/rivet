// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief q^2 in D0 -> pi- mu+ nu_mu and D+ -> pi0  mu+ nu_mu decays
  class BESIII_2018_I1655158 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1655158);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2_D0, 1, 1, 1);
      book(_h_q2_Dp, 2, 1, 1);
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
      // Loop over D mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::D0 or
									       Cuts::pid==PID::DPLUS )) {
        if (p.pid()==PID::D0 &&
	    isSemileptonicDecay(p, {PID::PIMINUS, PID::ANTIMUON, PID::NU_MU}) ) {
	  _h_q2_D0->fill(q2(p, PID::PIMINUS));
        }
	else if(p.pid()==PID::DPLUS &&
		isSemileptonicDecay(p, {PID::PI0, PID::ANTIMUON, PID::NU_MU}) ) {
	  _h_q2_Dp->fill(q2(p, PID::PI0));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_q2_D0, 1.);
      normalize(_h_q2_Dp, 1.);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_q2_D0, _h_q2_Dp;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1655158);


}
