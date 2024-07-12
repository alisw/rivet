// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief 
  class LHCB_2015_I1401225 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2015_I1401225);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      declare(D0, "D0");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, { 211,1}, {13,1}, {-13,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 321,1}, {-211,1}, {13,1}, {-13,1} };
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.decaying()[ix].pid()>0 && D0.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (D0.decaying()[ix].pid()<0 && D0.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
       	const Particle & Kp = D0.decayProducts()[ix].at(-sign*321)[0];
      	const Particle & pim= D0.decayProducts()[ix].at( sign*211)[0];
      	const Particle & mup= D0.decayProducts()[ix].at( 13)[0];
      	const Particle & mum= D0.decayProducts()[ix].at(-13)[0];
	_h[0]->fill((Kp.momentum()+pim.momentum()).mass()/MeV);
	_h[1]->fill((mum.momentum()+mup.momentum()).mass()/MeV);
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


  RIVET_DECLARE_PLUGIN(LHCB_2015_I1401225);

}
