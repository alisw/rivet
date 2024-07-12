// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Lambda Lambdabar eta
  class BESIII_2022_I2167804 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2167804);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      psi.addStable(PID::OMEGA);
      psi.addStable( 3122);
      psi.addStable(-3122);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,2.,7.,50,2.,7.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 3122,1}, {-3122,1}, { 221,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particle & eta = psi.decayProducts()[ix].at( 221)[0];
	const Particle & lam   = psi.decayProducts()[ix].at( 3122)[0];
	const Particle & lbar  = psi.decayProducts()[ix].at(-3122)[0];
	double mminus = (lbar.momentum()+eta.momentum()).mass2();
	double mplus  = (lam .momentum()+eta.momentum()).mass2();
	_h[0]->fill(sqrt(mplus ));
	_h[1]->fill(sqrt(mminus));
	_h[2]->fill((lam .momentum()+lbar.momentum()).mass());
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix]);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2167804);

}
