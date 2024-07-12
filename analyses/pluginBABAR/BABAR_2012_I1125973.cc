// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> pi, eta and omega decays
  class BABAR_2012_I1125973 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1125973);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==PID::BPLUS || Cuts::pid==PID::B0), "UFS");

      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_B0_pi[ix], 1, 1, ix+1);
	book(_h_Bp_pi[ix], 2, 1, ix+1);
	book(_nB[ix],"TMP/nB_"+toString(ix+1));
      }
      book(_h_omega,3,1,1);
      book(_h_eta  ,4,1,1);
    }

    // Calculate the Q2 using mother and daughter charged lepton
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
      // Get B+ Mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	if( p.pid()==PID::B0) {
	  _nB[0]->fill();
	  if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E}) ||
	      isSemileptonicDecay(p, {PID::PIMINUS, PID::ANTIMUON, PID::NU_MU})) {
	    double Q2 = q2(p,PID::PIMINUS);
            _h_B0_pi[0]->fill(Q2);
            _h_B0_pi[1]->fill(Q2);
	  }
	}
	else {
	  _nB[1]->fill();
	  if (isSemileptonicDecay(p, {PID::PI0, PID::POSITRON, PID::NU_E}) ||
	      isSemileptonicDecay(p, {PID::PI0, PID::ANTIMUON, PID::NU_MU})) {
	    double Q2 = q2(p,PID::PI0);
            _h_Bp_pi[0]->fill(Q2);
            _h_Bp_pi[1]->fill(Q2);
	  }
	  else if (isSemileptonicDecay(p, {PID::OMEGA, PID::POSITRON, PID::NU_E}) ||
		   isSemileptonicDecay(p, {PID::OMEGA, PID::ANTIMUON, PID::NU_MU})) {
	    _h_omega->fill(q2(p,PID::OMEGA));
	  }
	  else if (isSemileptonicDecay(p, {PID::ETA, PID::POSITRON, PID::NU_E}) ||
		   isSemileptonicDecay(p, {PID::ETA, PID::ANTIMUON, PID::NU_MU})) {
            _h_eta->fill(q2(p,PID::ETA));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // BR in units of 10^{-7}
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_B0_pi[ix], 1e7/ *_nB[0]);
	scale(_h_Bp_pi[ix], 1e7/ *_nB[1]);
      }
      scale(_h_omega,1e7/ *_nB[1]);
      scale(_h_eta  ,1e7/ *_nB[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nB[2];
    Histo1DPtr _h_B0_pi[2],_h_Bp_pi[2],_h_omega,_h_eta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1125973);

}
