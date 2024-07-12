// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> pi+pi-pi0
  class BES_2005_I689969 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_2005_I689969);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi2S(ufs);
      psi2S.addStable(PID::PI0);
      psi2S.addStable(PID::K0S);
      psi2S.addStable(PID::ETA);
      psi2S.addStable(PID::ETAPRIME);
      declare(psi2S, "psi2S");
      book(_h_pipi,1,1,1);
      book(_dalitz, "dalitz",50,0.,14.,50,0.0,14.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 111,1} };
      DecayedParticles psi2S = apply<DecayedParticles>(event, "psi2S");
      // loop over particles
      for(unsigned int ix=0;ix<psi2S.decaying().size();++ix) {
	if(!psi2S.modeMatches(ix,3,mode)) continue;
	const Particles & pi0 = psi2S.decayProducts()[ix].at( 111);
	const Particles & pip = psi2S.decayProducts()[ix].at( 211);
	const Particles & pim = psi2S.decayProducts()[ix].at(-211);
	double mminus = (pim[0].momentum()+pi0[0].momentum()).mass2();
	double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	double mneut  = (pip[0].momentum()+pim[0].momentum()).mass2();
	_h_pipi->fill(sqrt(mneut ));
	_h_pipi->fill(sqrt(mplus ));
	_h_pipi->fill(sqrt(mminus));
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_2005_I689969);

}
