// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> p Lambdabar0 pi- and Bs0 -> p Lambdabar0 K-
  class LHCB_2017_I1596893 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2017_I1596893);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==531);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 3122);
      B0.addStable(-3122);
      declare(B0, "B0");
      // book histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 2212,1},{-3122,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-2212,1},{ 3122,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 2212,1},{-3122,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-2212,1},{ 3122,1}, { 321,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign = 1, imode = 0;
	if(B0.decaying()[ix].abspid()==511) {
	  if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode1)) {
	    sign=1;
	  }
	  else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode1CC)) {
	    sign=-1;
	  }
	  else
	    continue;
	  imode=0;
	}
	else {
	  if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode2)) {
	    sign=1;
	  }
	  else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode2CC)) {
	    sign=-1;
	  }
	  else
	    continue;
	  imode=1;
	}
       	const Particle & pp     = B0.decayProducts()[ix].at( sign*2212)[0];
       	const Particle & LamBar = B0.decayProducts()[ix].at(-sign*3122)[0];
       	const Particle & pim    = B0.decayProducts()[ix].at(-sign*(imode==0 ? 211 : 321))[0];
       	_h[imode][0]->fill((pp .momentum()+LamBar.momentum()).mass());
       	_h[imode][1]->fill((pp .momentum()+pim   .momentum()).mass());
       	_h[imode][2]->fill((pim.momentum()+LamBar.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2017_I1596893);

}
