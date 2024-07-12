// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> D_s(*) K pi
  class BELLE_2009_I816632 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I816632);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 431);
      BP.addStable(-431);
      BP.addStable( 433);
      BP.addStable(-433);
      declare(BP, "BP");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_mass[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { -431,1}, { 321,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {  431,1}, {-321,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode2   = { { -433,1}, { 321,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {  433,1}, {-321,1}, {-211,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign=1, mode=0;
      	if      (BP.modeMatches(ix,3,mode1  )) {
	  mode = 0;
	  sign = 1;
	}
      	else if (BP.modeMatches(ix,3,mode1CC)) {
	  sign =-1;
	  mode = 0;
	}
	else if      (BP.modeMatches(ix,3,mode2  )) {
	  mode = 1;
	  sign = 1;
	}
      	else if (BP.modeMatches(ix,3,mode2CC)) {
	  sign =-1;
	  mode = 1;
	}
      	else continue;
	const Particle & Ds  = BP.decayProducts()[ix].at(-sign*(431+2*mode))[0];
       	const Particle & Kp  = BP.decayProducts()[ix].at( sign*321)[0];
       	_h_mass[mode]->fill((Ds.momentum()+Kp.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_mass[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I816632);

}
