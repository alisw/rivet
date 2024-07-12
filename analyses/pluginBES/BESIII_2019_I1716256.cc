// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> J/psi mu+mu-
  class BESIII_2019_I1716256 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1716256);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443||Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable(PID::JPSI);
      declare(chi, "chi");
      // Initialise and register projections
      for(unsigned int ix=0;ix<2;++ix) {
      	for(unsigned int iy=0;iy<2;++iy) {
      	  book(_h[ix][iy],1+ix,1,1+iy);
      	  book(_n[ix][iy],"TMP/n_"+toString(ix+1)+"_"+toString(iy+1));
      	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static double alpha= 7.2973525664e-3;
      static const map<PdgId,unsigned int> & mode1 = { { 443,1}, {22,1} };
      static const map<PdgId,unsigned int> & mode2 = { { 443,1}, {13,1}, {-13,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	unsigned int iloc=chi.decaying()[ix].pid()==20443 ? 0 : 1;
	_n[iloc][0]->fill();
	if(chi.modeMatches(ix,2,mode1)) {
	  _n[iloc][1]->fill();
	}
	else if(chi.modeMatches(ix,3,mode2)) {
	  const Particle & mum  = chi.decayProducts()[ix].at( 13)[0];
	  const Particle & mup  = chi.decayProducts()[ix].at(-13)[0];
	  const Particle & Jpsi = chi.decayProducts()[ix].at(443)[0];
	  double q = (mum.momentum()+mup.momentum()).mass();
	  // br
	  _h[iloc][0]->fill(q);
	  double M2 = chi.decaying()[ix].mass2();
	  double m  = Jpsi.mass();
	  // factor for form factor (eqn 2 arXiv:1709.02444)
	  double fact = alpha/3./M_PI/sqr(q)*
	    sqrt((1.-sqr(m+q)/M2)*(1.-sqr(m-q)/M2))/(1.-sqr(m)/M2)*
	    (1.+2.*sqr(mum.mass()/q))*sqrt(1.-4.*sqr(mum.mass()/q));
	  // additional factor of 2q from dq^2 = 2q dq
	  _h[iloc][1]->fill(q,1./fact/2/q);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h[ix][iy], 1./ *_n[ix][iy]);
	}
	scale(_h[ix][0], 1e5);
      }
    }
    /// @}

    /// @name Histograms
    ///@{
    Histo1DPtr _h[2][2];
    CounterPtr _n[2][2];

    ///@}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1716256);

}
