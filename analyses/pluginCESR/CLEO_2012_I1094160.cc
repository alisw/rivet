// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 K+/- pi-/+
  class CLEO_2012_I1094160 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2012_I1094160);


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
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_K0Km [ix],1,1,1+3*ix);
	book(_h_K0pip[ix],1,1,2+3*ix);
	book(_h_Kmpip[ix],1,1,3+3*ix);
	book(_h_K0Kp [ix],2,1,1+3*ix);
	book(_h_K0pim[ix],2,1,2+3*ix);
	book(_h_Kppim[ix],2,1,3+3*ix);
	book(_dalitz [ix],"dalitz_"+toString(ix+1),50,0.3,2.0,50,0.3,2.);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 310,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( !D0.modeMatches(ix,3,mode  ) &&
	    !D0.modeMatches(ix,3,modeCC) ) continue;
	const Particles & K0 = D0.decayProducts()[ix].at(310);
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & pip= D0.decayProducts()[ix].find( sign*211) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at( sign*211);
	const Particles & pim= D0.decayProducts()[ix].find(-sign*211) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at(-sign*211);
	const Particles & Kp = D0.decayProducts()[ix].find( sign*321) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at( sign*321);
	const Particles & Km = D0.decayProducts()[ix].find(-sign*321) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at(-sign*321);
       	// K0S K- pi+
       	if( Km.size()==1 && pip.size()==1) {
       	  double mK0pip = (K0[0].momentum()+pip[0].momentum() ).mass2();
       	  double mKmpip = (Km[0].momentum()+pip[0].momentum() ).mass2();
       	  double mKK    = (K0[0].momentum()+Km [0].momentum() ).mass2();
       	  for(unsigned int ix=0;ix<2;++ix) {
       	    _h_K0Km [ix]->fill(mKK   );
       	    _h_K0pip[ix]->fill(mK0pip);
       	    _h_Kmpip[ix]->fill(mKmpip);
       	  }
       	  _dalitz[0]->fill(mKmpip,mK0pip); 
       	}
       	// K0S K+ pi-
       	else if( Kp.size()==1 && pim.size()==1) {
       	  double mK0pim = (K0[0].momentum()+pim[0].momentum() ).mass2();
       	  double mKppim = (Kp[0].momentum()+pim[0].momentum() ).mass2();
       	  double mKK    = (K0[0].momentum()+Kp [0].momentum() ).mass2();
       	  for(unsigned int ix=0;ix<2;++ix) {
       	    _h_K0Kp [ix]->fill(mKK   );
       	    _h_K0pim[ix]->fill(mK0pim);
       	    _h_Kppim[ix]->fill(mKppim);
       	  }
       	  _dalitz[1]->fill(mKppim,mK0pim); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_K0Km [ix]);
	normalize(_h_K0pip[ix]);
	normalize(_h_Kmpip[ix]);
	normalize(_h_K0Kp [ix]);
	normalize(_h_K0pim[ix]);
	normalize(_h_Kppim[ix]);
	normalize(_dalitz [ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kmpip[2], _h_K0pip[2], _h_K0Km[2];
    Histo1DPtr _h_Kppim[2], _h_K0pim[2], _h_K0Kp[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2012_I1094160);

}
