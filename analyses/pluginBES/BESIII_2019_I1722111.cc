// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


  /// @brief psi(2S) -> p pbar phi
  class BESIII_2019_I1722111 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1722111);


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
      psi.addStable(PID::PHI);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,3.5,8,50,3.5,8);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 333,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particle & phi  = psi.decayProducts()[ix].at( 333)[0];
	const Particle & pp   = psi.decayProducts()[ix].at( 2212)[0];
	const Particle & pbar = psi.decayProducts()[ix].at(-2212)[0];
	if(phi.mass()<1.005 || phi.mass()>1.035) continue;
	double mminus = (pbar.momentum()+phi .momentum()).mass2();
	double mplus  = (pp  .momentum()+phi .momentum()).mass2();
	double mneut  = (pp  .momentum()+pbar.momentum()).mass2();
	_h[0]->fill(sqrt(mplus ));
	_h[1]->fill(sqrt(mplus ));
	_h[2]->fill(sqrt(mneut ));
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


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1722111);

}
