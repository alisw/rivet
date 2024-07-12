// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> phi phi K+
  class BELLE_2021_I1841899 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1841899);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::PI0);
      BP.addStable(PID::K0S);
      BP.addStable(PID::ETA);
      BP.addStable(PID::ETAPRIME);
      BP.addStable(PID::PHI);
      declare(BP, "BP");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{ 333,2}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 333,2}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	int sign = 1;
	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & Kp  = BP.decayProducts()[ix].at( sign*321);
	const Particles & phi = BP.decayProducts()[ix].at(      333);
	_h[0]->fill((phi[0].momentum()+phi[1].momentum()).mass());
	_h[1]->fill((Kp [0].momentum()+phi[0].momentum()).mass());
	_h[1]->fill((Kp [0].momentum()+phi[1].momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1841899);

}
