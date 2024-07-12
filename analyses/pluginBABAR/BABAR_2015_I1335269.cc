// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D D K
  class BABAR_2015_I1335269 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1335269);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 411);
      BB.addStable(-411);
      BB.addStable( 421);
      BB.addStable(-421);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { {-411,1},{ 421,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 411,1},{-421,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode2   = { {-421,1},{ 421,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 421,1},{-421,1}, {-321,1}};
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
	const Particle & D0  = BB.decayProducts()[ix].at( sign*421)[0];
	const Particle & Kp  = BB.decayProducts()[ix].at( sign*321)[0];
	int iD = imode==0 ? 411 : 421;
	const Particle & Dm  = BB.decayProducts()[ix].at(-sign*iD )[0];
	_h[imode][0]->fill((D0.momentum()+Dm.momentum()).mass());
	_h[imode][1]->fill((Kp.momentum()+Dm.momentum()).mass());
	_h[imode][2]->fill((Kp.momentum()+D0.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  normalize(_h[ix][iy], 1., false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1335269);

}
