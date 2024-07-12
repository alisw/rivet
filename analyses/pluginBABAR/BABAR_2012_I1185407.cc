// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> 5 and 5 pion final states
  class BABAR_2012_I1185407 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1185407);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      declare(TAU, "TAU");
      // histogram
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h_6pi[ix],3,1,1+ix);
	if(ix>2) continue;
	book(_h_3pi[ix],1,1,1+ix);
	book(_h_5pi[ix],2,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { {-211,2}, { 211,1},{111,3},{ 16,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 211,2}, {-211,1},{111,3},{-16,1}};
      static const map<PdgId,unsigned int> & mode2   = { {-211,3}, { 211,2},{ 16,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 211,3}, {-211,2},{-16,1}};
      static const map<PdgId,unsigned int> & mode3   = { {-211,3}, { 211,2},{111,1},{ 16,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 211,3}, {-211,2},{111,1},{-16,1}};
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for(unsigned int ix=0;ix<TAU.decaying().size();++ix) {
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
	int imode=-1;
      	if(TAU.modeMatches(ix,7,mode1  ) ||
	   TAU.modeMatches(ix,7,mode1CC))      imode = 0;
	else if(TAU.modeMatches(ix,6,mode2  ) ||
		TAU.modeMatches(ix,6,mode2CC)) imode = 1;
	else if(TAU.modeMatches(ix,7,mode3  ) ||
		TAU.modeMatches(ix,7,mode3CC)) imode = 2;
	else continue;
      	const Particles & pip = TAU.decayProducts()[ix].at( 211*sign);
	const Particles & pim = TAU.decayProducts()[ix].at(-211*sign);
	FourMomentum ptotal;
	for(unsigned int ix=0;ix<pim.size();++ix) ptotal +=pim[ix].momentum();
	for(unsigned int ix=0;ix<pip.size();++ix) ptotal +=pip[ix].momentum();
	if(imode==1) {
	  _h_5pi[2]->fill(ptotal.mass());
	  for(unsigned int ix=0;ix<pim.size();++ix) {
	    _h_5pi[1]->fill((ptotal-pim[ix].momentum()).mass());
	    for(unsigned int iy=0;iy<pip.size();++iy) {
	      _h_5pi[0]->fill((pim[ix].momentum()+pip[iy].momentum()).mass());
	    }
	  }
	}
	else {
	  const Particles  & pi0 = TAU.decayProducts()[ix].at(111);
	  FourMomentum ppi0;
	  for(unsigned int ix=0;ix<pi0.size();++ix) ppi0 +=pi0[ix].momentum();
	  ptotal+=ppi0;
	  if(imode==0 ) {
	    _h_3pi[0]->fill(ppi0.mass());
	    _h_3pi[2]->fill(ptotal.mass());
	    for(unsigned int ix=0;ix<pim.size();++ix) {
	      for(unsigned int iy=0;iy<pi0.size();++iy) {
		_h_3pi[1]->fill((pip[0].momentum()+pim[ix].momentum()+pi0[iy].momentum()).mass());
	      }
	    }	    
	  }
	  else {
	    _h_6pi[3]->fill(ptotal.mass());
	    for(unsigned int ix=0;ix<pim.size();++ix) {
	      _h_6pi[2]->fill((ptotal-pim[ix].momentum()).mass());
	      for(unsigned int iy=0;iy<pip.size();++iy) {
		_h_6pi[0]->fill((pim[ix].momentum()+pip[iy].momentum()).mass());
		_h_6pi[1]->fill((pim[ix].momentum()+pip[iy].momentum()+pi0[0].momentum()).mass());
	      }
	    }
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix) {
	normalize(_h_6pi[ix]);
	if(ix>2) continue;
	normalize(_h_3pi[ix]);
	normalize(_h_5pi[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_3pi[3],_h_5pi[3],_h_6pi[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1185407);

}
