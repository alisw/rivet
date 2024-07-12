// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  J/psi -> phi eta eta'
  class BESIII_2019_I1712729 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1712729);


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
      psi.addStable(PID::OMEGA);
      psi.addStable(PID::PHI);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,3.5,7.,50,2.0,5.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 333,1}, { 221,1}, { 331,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particle & phi  = psi.decayProducts()[ix].at( 333)[0];
	const Particle & eta  = psi.decayProducts()[ix].at( 221)[0];
	const Particle & etap = psi.decayProducts()[ix].at( 331)[0];
	double mphieta  = (phi.momentum()+eta .momentum()).mass2();
	double mphietap = (phi.momentum()+etap.momentum()).mass2();
	_h[0]->fill(sqrt(mphietap));
	_h[1]->fill(sqrt(mphietap));
	_dalitz->fill(mphietap,mphieta);
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


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1712729);

}
