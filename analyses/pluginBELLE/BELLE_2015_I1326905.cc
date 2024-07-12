// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 ->\to Ds- KS0 pi+ and B+ -> Ds- K+ K^+ 
  class BELLE_2015_I1326905 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1326905);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 431);
      BB.addStable(-431);
      BB.addStable( 310);
      declare(BB, "BB");
      //histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { {-431,1},{ 310,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 431,1},{ 310,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode2   = { {-431,1},{ 321,2}};
      static const map<PdgId,unsigned int> & mode2CC = { { 431,1},{-321,2}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int sign = 1, imode = 0;
      	if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	  imode=0;
      	  sign=1;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
	  imode=0;
      	  sign=-1;
      	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	  imode=1;
      	  sign=1;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
	  imode=1;
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & Ds  = BB.decayProducts()[ix].at(-sign*431)[0];
	FourMomentum pK;
	if(imode==0) {
	  pK = BB.decayProducts()[ix].at(310)[0].momentum();
	}
	else {
	  const Particles & Kp = BB.decayProducts()[ix].at(sign*321);
	  pK = Kp[0].momentum().p3().mod()>Kp[1].momentum().p3().mod() ?
	    Kp[1].momentum() : Kp[0].momentum();
	}
	_h[imode]->fill((Ds.momentum()+pK).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h[ix], 1., false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1326905);

}
