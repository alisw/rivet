// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+ pi0 pi0
  class BESIII_2022_I2088337 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2088337);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");
      // histograms
      for(unsigned int ix=0;ix<7;++ix) {
	book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {310,1}, { 211,1}, { 111,2} };
      static const map<PdgId,unsigned int> & modeCC = { {310,1}, {-211,1}, { 111,2} };
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
	if ( !DD.modeMatches(ix,4,mode) && !DD.modeMatches(ix,4,modeCC)) continue;
	const Particles & KS0 = DD.decayProducts()[ix].at(      310);
	const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	_h[0]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	_h[3]->fill((KS0[0].momentum()+pip[0].momentum()).mass());
	_h[5]->fill((KS0[0].momentum()+pi0[0].momentum()+pi0[1].momentum()).mass());
	_h[6]->fill((pip[0].momentum()+pi0[0].momentum()+pi0[1].momentum()).mass());
	for(unsigned int ix=0;ix<2;++ix) {
	  _h[1]->fill((KS0[0].momentum()+pi0[ix].momentum()).mass());
	  _h[2]->fill((pip[0].momentum()+pi0[ix].momentum()).mass());
	  _h[4]->fill((KS0[0].momentum()+pip[0].momentum()+pi0[ix].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<7;++ix) {
	normalize(_h[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2088337);

}
