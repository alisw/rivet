// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> K+ pi+ pi- gamma
  class BABAR_2016_I1409292 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2016_I1409292);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::PI0);
      declare(BP, "BP");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix],1+ix,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{ 211,1}, {-211,1}, {22,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, {-211,1}, {22,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,4,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,4,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & Kp  = BP.decayProducts()[ix].at( sign*321)[0];
	const Particle & pip = BP.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim = BP.decayProducts()[ix].at(-sign*211)[0];
	_h[0]->fill((Kp.momentum()+pim.momentum()+pip.momentum()).mass());
	_h[1]->fill((Kp.momentum()+pim.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2016_I1409292);

}
