// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> Lambda_c+ Lambdabar_c- K-
  class BELLE_2018_I1642436 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1642436);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 4122);
      BP.addStable(-4122);
      declare(BP, "BP");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4122,1},{-4122,1}, {-321,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 4122,1},{-4122,1}, { 321,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
       	const Particle & LamC    = BP.decayProducts()[ix].at( sign*4122)[0];
       	const Particle & LamCBar = BP.decayProducts()[ix].at(-sign*4122)[0];
       	const Particle & Km      = BP.decayProducts()[ix].at(-sign*321 )[0];
       	_h[0]->fill((LamC.momentum()+Km     .momentum()).mass());
      	_h[1]->fill((LamC.momentum()+LamCBar.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1642436);

}
