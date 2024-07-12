// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta -> pi0 gamma gamma
  class A2_2014_I1297221 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(A2_2014_I1297221);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==PID::ETA);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      declare(ETA, "ETA");
      // histos
      for(unsigned int ix=0;ix<3;++ix) book(_h[ix],1,1,1+ix);
      book(_c,"TMP/nEta");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {111,1}, { 22,2} };
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	_c->fill();
	// select right decay mode
	if ( !ETA.modeMatches(ix,3,mode)) continue;
	const Particles & gam = ETA.decayProducts()[ix].at(22);
	double mass2 = (gam[0].momentum()+gam[1].momentum()).mass2();
	for(unsigned int ix=0;ix<3;++ix) _h[ix]->fill(mass2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // eta width in eV from PDG2022
      double gam = 1.31e3;
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h[ix],gam / *_c);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(A2_2014_I1297221);

}
