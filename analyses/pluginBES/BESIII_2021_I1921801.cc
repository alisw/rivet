// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> omega KS0 KS0
  class BESIII_2021_I1921801 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1921801);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi2S(ufs);
      psi2S.addStable(PID::PI0);
      psi2S.addStable(PID::K0S);
      psi2S.addStable(PID::ETA);
      psi2S.addStable(PID::ETAPRIME);
      psi2S.addStable(PID::OMEGA);
      declare(psi2S, "psi2S");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz" ,50,1.,11.,50,1.,11.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { { 310,2},{ 223,1}};
      DecayedParticles psi2S = apply<DecayedParticles>(event, "psi2S");
      // loop over particles
      for(unsigned int ix=0;ix<psi2S.decaying().size();++ix) {
	if(!psi2S.modeMatches(ix,3,mode)) continue;
	const Particles & K0    = psi2S.decayProducts()[ix].at(310);
	const Particles & omega = psi2S.decayProducts()[ix].at(223);
	double m2Komega[2]={(K0[0].momentum()+omega[0].momentum()).mass2(),
			    (K0[1].momentum()+omega[0].momentum()).mass2()};
	for(unsigned int ix=0;ix<2;++ix) {
	  double mKo = sqrt(m2Komega[ix]);
	  _h[0]->fill(mKo);
	  _h[1]->fill(mKo);
	}
	_h[2]->fill((K0[0].momentum()+K0[1].momentum()).mass());
	_dalitz->fill(m2Komega[0],m2Komega[1]);
	_dalitz->fill(m2Komega[1],m2Komega[0]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1921801);

}
