// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> eta', eta, pi0 e+e-
  class BESIII_2014_I1287631 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2014_I1287631);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==PID::JPSI);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "PSI");
      for(unsigned int ix=0;ix<5;++ix)
	book(_h[ix], 1, 1, ix+1);
      book(_h[5],2,1,1);
      book(_c,"TMP/den");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode0   = { {PID::ETAPRIME,1},{ PID::PHOTON,1} };
      static const map<PdgId,unsigned int> & mode1   = { {PID::ETAPRIME,1},{ 11,1}, { -11,1}};
      static const map<PdgId,unsigned int> & mode2   = { {PID::ETA     ,1},{ 11,1}, { -11,1}};
      static const map<PdgId,unsigned int> & mode3   = { {PID::PI0     ,1},{ 11,1}, { -11,1}};
      DecayedParticles psi = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	unsigned int imode=0;
	if(psi.modeMatches(ix,2,mode0)) {
	  _c->fill();
	  continue;
	}
	else if(psi.modeMatches(ix,3,mode1)) {
	  imode=1;
	}
	else if(psi.modeMatches(ix,3,mode2)) {
	  imode=2;
	}
	else if(psi.modeMatches(ix,3,mode3)) {
	  imode=3;
	}
	else
	  continue;
	// e+ e-
	const Particle & em = psi.decayProducts()[ix].at( 11)[0];
	const Particle & ep = psi.decayProducts()[ix].at(-11)[0];
	double q = (ep.momentum()+em.momentum()).mass();
	if(imode==1) {
	  _h[0]->fill(q);
	  _h[1]->fill(q);
	  double me = em.mass();
	  double beta = sqrt(1.-4.*sqr(me/q));
	  double mJpsi = psi.decaying()[ix].mass();
	  double mEta  = psi.decayProducts()[ix].at(331)[0].mass();
	  double p = sqrt(sqr(1.+sqr(q)/(sqr(mJpsi)-sqr(mEta)))-4.*sqr(mJpsi*q/(sqr(mJpsi)-sqr(mEta))));
	  double fact = beta*GeV/q*(1.+2.*sqr(me/q))*pow(p,3);
	  _h[5]->fill(q,1./fact);
	}
	else if(imode==2) {
	  _h[2]->fill(q);
	  _h[3]->fill(q);
	}
	else {
	  _h[4]->fill(q);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix],1.,false);
      static double alpha= 7.2973525664e-3;
      scale(_h[5], 1.5 *M_PI/alpha/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2014_I1287631);

}
