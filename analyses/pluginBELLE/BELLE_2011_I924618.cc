// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> pbar Lambda D0
  class BELLE_2011_I924618 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2011_I924618);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 3122);
      BP.addStable(-3122);
      BP.addStable( 421);
      BP.addStable(-421);
      declare(BP, "BP");
      // histos
      book(_h,1,1,1);
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-2212,1},{ 3122,1}, { 421,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 2212,1},{-3122,1}, {-421,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	_c->fill();
      	int sign = 1;
      	if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & pbar = BP.decayProducts()[ix].at(-sign*2212)[0];
	const Particle & Lam  = BP.decayProducts()[ix].at( sign*3122)[0];
	_h->fill((pbar.momentum()+Lam.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1e6/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2011_I924618);

}
