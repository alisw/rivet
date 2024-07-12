// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> p nbar pi- p nbar pi- pi0
  class BESIII_2012_I1128258 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1128258);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441 ||
						Cuts::pid==445   ||
						Cuts::pid==20444);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      declare(chi, "chi");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_h[3],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 2212,1}, {-2112,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-2212,1}, { 2112,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 2212,1}, {-2112,1}, {-211,1}, {111,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-2212,1}, { 2112,1}, { 211,1}, {111,1}};
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
      	int mode=-1,sign=0;
      	if (chi.decaying()[ix].pid()==10441) {
      	  if(chi.modeMatches(ix,3,mode1)) {
      	    mode=0;
      	    sign=1;
      	  }
      	  else if(chi.modeMatches(ix,3,mode1CC)) {
      	    mode= 0;
      	    sign=-1;
      	  }
	

       	}
      	if(mode==0) {
      	  const Particle & pim  = chi.decayProducts()[ix].at( -211*sign)[0];
      	  const Particle & pp   = chi.decayProducts()[ix].at( 2212*sign)[0];
      	  const Particle & nbar = chi.decayProducts()[ix].at(-2112*sign)[0];
      	  _h[0]->fill((pp  .momentum()+pim .momentum()).mass());
      	  _h[1]->fill((nbar.momentum()+pim .momentum()).mass());
      	  _h[2]->fill((pp  .momentum()+nbar.momentum()).mass());
	  continue;
      	}
      	else if(chi.modeMatches(ix,4,mode2)) {
      	  mode=1;
      	  sign=1;
      	}
      	else if(chi.modeMatches(ix,4,mode2CC)) {
      	  mode= 1;
      	  sign=-1;
      	}
      	else
      	  continue;
      	const Particle & pim = chi.decayProducts()[ix].at( -211*sign)[0];
      	const Particle & pi0 = chi.decayProducts()[ix].at(  111     )[0];
      	_h[3]->fill((pim.momentum()+pi0.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1128258);

}
