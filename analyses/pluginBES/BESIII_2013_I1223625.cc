// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> p pbar omega
  class BESIII_2013_I1223625 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1223625);


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
      declare(psi, "psi");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,2.5,5.,50,2.5,5.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 223,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & omega = psi.decayProducts()[ix].at( 223);
	const Particles & pp    = psi.decayProducts()[ix].at( 2212);
	const Particles & pbar  = psi.decayProducts()[ix].at(-2212);
	double mminus = (pbar[0].momentum()+omega[0].momentum()).mass2();
	double mplus  = (pp  [0].momentum()+omega[0].momentum()).mass2();
	double mneut  = (pp  [0].momentum()+pbar [0].momentum()).mass2();
	_h[0]->fill(sqrt(mplus ));
	_h[1]->fill(sqrt(mneut ));
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


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1223625);

}
