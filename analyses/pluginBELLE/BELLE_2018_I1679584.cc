// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> Lambda_c+ Lambdabar_c- Kbar0
  class BELLE_2018_I1679584 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1679584);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 4122);
      B0.addStable(-4122);
      B0.addStable( 310);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4122,1},{-4122,1}, { 310,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	if(!B0.modeMatches(ix,3,mode))
      	  continue;
	const Particle & LamC    = B0.decayProducts()[ix].at( 4122)[0];
	const Particle & LamCBar = B0.decayProducts()[ix].at(-4122)[0];
	const Particle & K0      = B0.decayProducts()[ix].at(  310)[0];
	if(B0.decaying()[ix].pid()<0)
	  _h[0]->fill((LamC   .momentum()+K0.momentum()).mass());
	else
	  _h[0]->fill((LamCBar.momentum()+K0.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1679584);

}
