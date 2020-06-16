// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D -> eta semi-leptonic q^2
  class CLEOC_2011_I875526 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2011_I875526);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_q2_Dp_eta_A,1,1,1);
      book(_h_q2_Dp_eta_B,1,1,2);
      book(_nDp,"TMP/nDp");
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
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::DPLUS )) {
	_nDp->fill();
	if(isSemileptonicDecay(p, {PID::ETA, PID::POSITRON, PID::NU_E})  ||
	   isSemileptonicDecay(p, {PID::ETA, PID::ELECTRON, PID::NU_EBAR})) {
	  double q = q2(p, PID::ETA);
	  _h_q2_Dp_eta_A->fill(q);
	  _h_q2_Dp_eta_B->fill(q);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_q2_Dp_eta_A, 1./ *_nDp);
      scale(_h_q2_Dp_eta_B, 1./ *_nDp);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_q2_Dp_eta_A,_h_q2_Dp_eta_B;
    CounterPtr _nD0,_nDp;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOC_2011_I875526);

}
