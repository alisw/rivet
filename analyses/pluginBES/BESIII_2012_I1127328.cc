// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> K+ K- pi0 
  class BESIII_2012_I1127328 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1127328);


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
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,0.,11.,50,0.0,11.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, {-321,1}, { 111,1} };
      DecayedParticles psi2S = apply<DecayedParticles>(event, "psi2S");
      // loop over particles
      for(unsigned int ix=0;ix<psi2S.decaying().size();++ix) {
	if(!psi2S.modeMatches(ix,3,mode)) continue;
	const Particles & pi0 = psi2S.decayProducts()[ix].at( 111);
	const Particles & Kp  = psi2S.decayProducts()[ix].at( 321);
	const Particles & Km  = psi2S.decayProducts()[ix].at(-321);
	double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	double mneut  = (Kp[0].momentum()+Km [0].momentum()).mass2();
	_h[0]->fill(sqrt(mneut ));
	_h[1]->fill(sqrt(mplus ));
	_h[1]->fill(sqrt(mminus));
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


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1127328);

}
