// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> Sigma_c++ pbar pi- pi-
  class BABAR_2012_I1127599 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1127599);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 4222);
      BP.addStable(-4222);
      declare(BP, "BP");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4222,1},{-2212,1}, {-211,2}};
      static const map<PdgId,unsigned int> & modeCC = { {-4222,1},{ 2212,1}, { 211,2}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,4,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,4,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle  & SigC = BP.decayProducts()[ix].at( sign*4222)[0];
	const Particle  & pbar = BP.decayProducts()[ix].at(-sign*2212)[0];
	const Particles & pim  = BP.decayProducts()[ix].at(-sign*211 );
	for(unsigned int ix=0;ix<2;++ix) {
	  _h[0]->fill((pbar.momentum()+pim[ix].momentum()).mass());
	  _h[1]->fill((SigC.momentum()+pim[ix].momentum()).mass());
	}
	_h[2]->fill((SigC.momentum()+pim[0].momentum()+pim[1].momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1127599);

}
