// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_b0 -> Lambda_c+ semi-leptonic
  class LHCB_2017_I1621811 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2017_I1621811);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2 , 1, 1, 1);
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
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==5122)) {
	if(p.pid()>0) {
	  if(isSemileptonicDecay(p, {4122, PID::MUON, PID::NU_MUBAR})) {
	    _h_q2->fill(q2(p, 4122));
	  }
	}
	else {
	  if(isSemileptonicDecay(p, {-4122, PID::ANTIMUON, PID::NU_MU})) {
	    _h_q2->fill(q2(p, -4122));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize to unity
      normalize(_h_q2);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_q2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2017_I1621811);

}
