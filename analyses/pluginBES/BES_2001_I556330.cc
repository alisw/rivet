// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> p pbar eta
  class BES_2001_I556330 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_2001_I556330);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
      book(_dalitz, "dalitz",50,2.,5.,50,2.,5.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 221,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & eta  = psi.decayProducts()[ix].at( 221);
	const Particles & pp   = psi.decayProducts()[ix].at( 2212);
	const Particles & pbar = psi.decayProducts()[ix].at(-2212);
	double mminus = (pbar[0].momentum()+eta[0].momentum()).mass2();
	double mplus  = (pp  [0].momentum()+eta[0].momentum()).mass2();
	_h[0]->fill(sqrt(mplus ));
	_h[0]->fill(sqrt(mminus));
	_h[1]->fill((pp  [0].momentum()+pbar[0].momentum()).mass());
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_dalitz);
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_2001_I556330);

}
