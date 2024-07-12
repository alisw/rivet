// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 ->  2pi+2pi-
  class FOCUS_2007_I741543 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2007_I741543);


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
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix   ],1,1,1+ix);
      book(_h[4],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode1   = { { 211,2}, { -211,2}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = D0.decaying()[ix].pid()/421;
	if ( D0.modeMatches(ix,4,mode1)) {
	  const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	  const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	  set<double> mpm; 
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      double m = (pip[ix].momentum()+pim[iy].momentum()).mass();
	      _h[0]->fill(m);
	      mpm.insert(m);
	    }
	  }
	  _h[1]->fill(*mpm.rbegin());
	  _h[2]->fill(*mpm.begin());
	  FourMomentum ppp = pip[0].momentum()+pip[1].momentum();
	  _h[3]->fill(ppp.mass());
	  FourMomentum pmm = pim[0].momentum()+pim[1].momentum();
	  _h[3]->fill(pmm.mass());
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h[4]->fill((ppp+pim[ix].momentum()).mass());
	    _h[4]->fill((pmm+pip[ix].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2007_I741543);

}
