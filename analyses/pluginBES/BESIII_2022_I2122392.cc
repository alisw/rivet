// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Lambda Lambdabar omega
  class BESIII_2022_I2122392 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2122392);


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
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,3.,7.,50,3.,7.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 3122,1}, {-3122,1}, { 223,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & omega = psi.decayProducts()[ix].at( 223);
	double mOmega=omega[0].mass();
	if(mOmega>.873 || mOmega<.753) continue;
	const Particles & lam   = psi.decayProducts()[ix].at( 3122);
	const Particles & lbar  = psi.decayProducts()[ix].at(-3122);
	double mminus = (lbar[0].momentum()+omega[0].momentum()).mass2();
	double mplus  = (lam [0].momentum()+omega[0].momentum()).mass2();
	_h[0]->fill(mplus );
	_h[1]->fill(mminus);
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix]);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2122392);

}
