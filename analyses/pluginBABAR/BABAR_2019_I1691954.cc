// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K- pi+ e+e-
  class BABAR_2019_I1691954 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2019_I1691954);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      declare(D0, "D0");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
      book(_h[2],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-321,1}, { 211,1}, {11,1}, {-11,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 321,1}, {-211,1}, {11,1}, {-11,1} };
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.decaying()[ix].pid()>0 && D0.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (D0.decaying()[ix].pid()<0 && D0.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
       	const Particle & Kp  = D0.decayProducts()[ix].at(-sign*321)[0];
      	const Particle & pim = D0.decayProducts()[ix].at( sign*211)[0];
      	const Particle & ep  = D0.decayProducts()[ix].at( 11)[0];
      	const Particle & em  = D0.decayProducts()[ix].at(-11)[0];
	double mee = (em.momentum()+ep.momentum()).mass();
	if(mee>0.675 && mee<0.875) {
	  _h[1]->fill((Kp.momentum()+pim.momentum()).mass());
	  _h[0]->fill(mee);
	}
	if(mee>.2) _h[2]->fill(mee);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2019_I1691954);

}
