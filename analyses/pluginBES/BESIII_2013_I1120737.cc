// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> p pbar pi0
  class BESIII_2013_I1120737 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1120737);


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
      declare(psi, "psi");
      book(_h_ppi0,1,1,1);
      book(_dalitz, "dalitz",50,1.,8.,50,1.,8.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 111,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & pi0  = psi.decayProducts()[ix].at( 111);
	const Particles & pp   = psi.decayProducts()[ix].at( 2212);
	const Particles & pbar = psi.decayProducts()[ix].at(-2212);
	double mminus = (pbar[0].momentum()+pi0 [0].momentum()).mass2();
	double mplus  = (pp  [0].momentum()+pi0 [0].momentum()).mass2();
	_h_ppi0->fill(sqrt(mplus ));
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_dalitz);
      normalize(_h_ppi0);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_ppi0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1120737);

}
