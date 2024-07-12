// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief BbarS0 -> Lambda_c+ Lambdabar0 pi-
  class BELLE_2013_I1230342 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1230342);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==531);
      declare(ufs, "UFS");
      DecayedParticles BS0(ufs);
      BS0.addStable( 3122);
      BS0.addStable(-3122);
      BS0.addStable( 4122);
      BS0.addStable(-4122);
      declare(BS0, "BS0");
      // histos
      book(_c,"TMP/nB");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4122,1},{-3122,1}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-4122,1},{ 3122,1}, { 211,1}};
      DecayedParticles BS0 = apply<DecayedParticles>(event, "BS0");
      // loop over particles
      for(unsigned int ix=0;ix<BS0.decaying().size();++ix) {
	_c->fill();
      	int sign = 1;
      	if (BS0.decaying()[ix].pid()<0 && BS0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (BS0.decaying()[ix].pid()>0 && BS0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & LamC   = BS0.decayProducts()[ix].at( sign*4122)[0];
	const Particle & LamBar = BS0.decayProducts()[ix].at(-sign*3122)[0];
	const Particle & pim    = BS0.decayProducts()[ix].at(-sign*211 )[0];
	_h[0]->fill((LamC.momentum()+LamBar.momentum()).mass());
	_h[1]->fill((LamC.momentum()+pim   .momentum()).mass());
	_h[2]->fill((pim .momentum()+LamBar.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	scale(_h[ix],1e4/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1230342);

}
