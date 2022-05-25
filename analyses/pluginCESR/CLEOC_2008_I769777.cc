// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D -> pi,K semi-leptonic q^2
  class CLEOC_2008_I769777 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2008_I769777);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_q2_D0_pi,1,1,1);
      book(_h_q2_Dp_pi,1,1,2);
      book(_h_q2_D0_K ,1,1,3);
      book(_h_q2_Dp_K ,1,1,4);
      book(_nD0,"TMP/nD0");
      book(_nDp,"TMP/nDp");
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::abspid==abs(mesonID))[0];
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
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::D0 or
									       Cuts::abspid==PID::DPLUS )) {
        if (p.abspid()==PID::D0) {
	  _nD0->fill();
	  if(isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E}) ||
	     isSemileptonicDecay(p, {PID::PIPLUS , PID::ELECTRON, PID::NU_EBAR}) )
	    _h_q2_D0_pi->fill(q2(p, PID::PIMINUS));
	  else if(isSemileptonicDecay(p, {PID::KMINUS, PID::POSITRON, PID::NU_E}) ||
		  isSemileptonicDecay(p, {PID::KPLUS , PID::ELECTRON, PID::NU_EBAR}))
	    _h_q2_D0_K ->fill(q2(p, PID::KMINUS));
        }
	else if(p.abspid()==PID::DPLUS) {
	  _nDp->fill();
	  if(isSemileptonicDecay(p, {PID::PI0, PID::POSITRON, PID::NU_E})  ||
	     isSemileptonicDecay(p, {PID::PI0, PID::ELECTRON, PID::NU_EBAR}))
	    _h_q2_Dp_pi->fill(q2(p, PID::PI0));
	  else if(isSemileptonicDecay(p, {-311, PID::POSITRON, PID::NU_E}))
	    _h_q2_Dp_K ->fill(q2(p, -311));
	  else if(isSemileptonicDecay(p, { 311, PID::ELECTRON, PID::NU_EBAR}))
	    _h_q2_Dp_K ->fill(q2(p, 311));
	  else if(isSemileptonicDecay(p, {PID::K0S, PID::POSITRON, PID::NU_E}) ||
		  isSemileptonicDecay(p, {PID::K0S, PID::ELECTRON, PID::NU_EBAR}))
	    _h_q2_Dp_K ->fill(q2(p, PID::K0S));
	  else if(isSemileptonicDecay(p, {PID::K0L, PID::POSITRON, PID::NU_E}) ||
		  isSemileptonicDecay(p, {PID::K0L, PID::ELECTRON, PID::NU_EBAR}))
	    _h_q2_Dp_K ->fill(q2(p, PID::K0L));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_q2_D0_pi,100./ *_nD0);
      scale(_h_q2_D0_K ,100./ *_nD0);
      scale(_h_q2_Dp_pi,100./ *_nDp);
      scale(_h_q2_Dp_K ,100./ *_nDp);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_q2_D0_pi, _h_q2_D0_K, _h_q2_Dp_pi, _h_q2_Dp_K;
    CounterPtr _nD0,_nDp;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I769777);

}
