// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_cJ -> Lambda Lambdabar0 eta
  class BESIII_2022_I2166668 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2166668);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441 or
						Cuts::pid==20443 or
						Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::ETA);
      chi.addStable( PID::LAMBDA);
      chi.addStable(-PID::LAMBDA);
      declare(chi, "chi");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 3122,1}, {-3122,1}, { 221,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,3,mode)) continue;
	const Particle & lam = chi.decayProducts()[ix].at( 3122)[0];
	const Particle & lamb = chi.decayProducts()[ix].at(-3122)[0];
	const Particle & eta  = chi.decayProducts()[ix].at(  221)[0];
	_h[0]->fill((lam .momentum()+lamb.momentum()).mass());
	_h[1]->fill((lam .momentum()+eta .momentum()).mass());
	_h[2]->fill((lamb.momentum()+eta .momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2166668);

}
