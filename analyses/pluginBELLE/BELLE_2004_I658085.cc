// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D* + pions
  class BELLE_2004_I658085 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I658085);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PI0);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 423);
      BB.addStable(-423);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      for (unsigned int ix=0;ix<6;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // decay modes
      static const map<PdgId,unsigned int> & mode1   = { { -413,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {  413,1}, {-211,2}, { 211,1} };
      static const map<PdgId,unsigned int> & mode2   = { { -413,1}, { 211,3}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2CC = { {  413,1}, {-211,3}, { 211,1} };
      static const map<PdgId,unsigned int> & mode3   = { { -413,1}, { 211,3}, {-211,2} };
      static const map<PdgId,unsigned int> & mode3CC = { {  413,1}, {-211,3}, { 211,2} };
      static const map<PdgId,unsigned int> & mode4   = { { -423,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode4CC = { {  423,1}, {-211,2}, { 211,1} };
      static const map<PdgId,unsigned int> & mode5   = { { -423,1}, { 211,2}, {-211,2} };
      static const map<PdgId,unsigned int> & mode5CC = { {  423,1}, {-211,2}, { 211,2} };
      static const map<PdgId,unsigned int> & mode6   = { { -423,1}, { 211,3}, {-211,2} };
      static const map<PdgId,unsigned int> & mode6CC = { {  423,1}, {-211,3}, { 211,2} };
      // loop over particles
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      int imode = -1;
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
	if ( (sign== 1 && BB.modeMatches(ix,4,mode1) ) ||
	     (sign==-1 && BB.modeMatches(ix,4,mode1CC) ) )
	  imode=0;
	else if ( (sign== 1 && BB.modeMatches(ix,5,mode2) ) ||
		  (sign==-1 && BB.modeMatches(ix,5,mode2CC) ) )
	  imode=1;
	else if ( (sign== 1 && BB.modeMatches(ix,6,mode3) ) ||
		  (sign==-1 && BB.modeMatches(ix,6,mode3CC) ) )
	  imode=2;
	else if ( (sign== 1 && BB.modeMatches(ix,4,mode4) ) ||
		  (sign==-1 && BB.modeMatches(ix,4,mode4CC) ) )
	  imode=3;
	else if ( (sign== 1 && BB.modeMatches(ix,5,mode5) ) ||
		  (sign==-1 && BB.modeMatches(ix,5,mode5CC) ) )
	  imode=4;
	else if ( (sign== 1 && BB.modeMatches(ix,6,mode6) ) ||
		  (sign==-1 && BB.modeMatches(ix,6,mode6CC) ) )
	  imode=5;
	else
	  continue;
	FourMomentum ptotal;
	for(const Particle & p : BB.decayProducts()[ix].at( sign*211) ) {
	  ptotal+=p.momentum();
	}
	for(const Particle & p : BB.decayProducts()[ix].at(-sign*211) ) {
	  ptotal+=p.momentum();
	}
	_h[imode]->fill(ptotal.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<6;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I658085);

}
