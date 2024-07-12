// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> K+K- eta
  class BESIII_2020_I1771616 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1771616);


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
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,0.,11.,50,0.0,11.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, {-321,1}, { 221,1} };
      DecayedParticles psi2S = apply<DecayedParticles>(event, "psi2S");
      // loop over particles
      for(unsigned int ix=0;ix<psi2S.decaying().size();++ix) {
	if(!psi2S.modeMatches(ix,3,mode)) continue;
	const Particles & eta = psi2S.decayProducts()[ix].at( 221);
	const Particles & Kp  = psi2S.decayProducts()[ix].at( 321);
	const Particles & Km  = psi2S.decayProducts()[ix].at(-321);
	double mminus = (Km[0].momentum()+eta[0].momentum()).mass2();
	double mplus  = (Kp[0].momentum()+eta[0].momentum()).mass2();
	double mneut  = (Kp[0].momentum()+Km [0].momentum()).mass2();
	_h[0]->fill(sqrt(mneut ));
	_h[1]->fill(sqrt(mplus ));
	_h[2]->fill(sqrt(mminus));
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


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1771616);

}
