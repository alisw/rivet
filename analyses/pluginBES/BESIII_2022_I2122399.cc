// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> p K- e+ nu_e
  class BESIII_2022_I2122399 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2122399);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==4122);
      declare(ufs, "UFS");
      DecayedParticles LambdaC(ufs);
      LambdaC.addStable(PID::PI0);
      LambdaC.addStable(PID::K0S);
      LambdaC.addStable(PID::ETA);
      LambdaC.addStable(PID::ETAPRIME);
      declare(LambdaC, "LambdaC");
      
      // Book histogram
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { -321,1}, { 2212,1}, {-11,1}, { 12,1}};
      DecayedParticles LambdaC = apply<DecayedParticles>(event, "LambdaC");
      // loop over particles
      for(unsigned int ix=0;ix<LambdaC.decaying().size();++ix) {
	if ( !LambdaC.modeMatches(ix,4,mode) ) continue;
       	const Particle & Km = LambdaC.decayProducts()[ix].at(-321 )[0];
       	const Particle & pp = LambdaC.decayProducts()[ix].at( 2212)[0];
	_h->fill((Km.momentum()+pp.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2122399);

}
