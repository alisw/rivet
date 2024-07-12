// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta' -> gamma gamma pi0
  class BESIII_2016_I1504943 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1504943);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==331);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      declare(ETA, "ETA");
      // Book histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_br[ix],1,1,1+ix);
      book(_h_m, 2, 1, 1);
      book(_netap, "TMP/netap");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {111,1}, { 22,2} };
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	_netap->fill();
	// select right decay mode
	if ( !ETA.modeMatches(ix,3,mode)) continue;
	const Particles & gam = ETA.decayProducts()[ix].at(22);
	double mass2 = (gam[0].momentum()+gam[1].momentum()).mass2();
	_h_m->fill(mass2);
	_h_br[0]->fill(.5);
	if(any(ETA.decaying()[ix].children(), hasAbsPID(PID::OMEGA))) _h_br[1]->fill(.5);
	if(ETA.decaying()[ix].children().size()==3) _h_br[2]->fill(0.5);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // eta' width in kev
      double gammaEtap = 0.188e3;
      scale(_h_m, gammaEtap/ *_netap);
      for(unsigned int ix=0;ix<3;++ix)
	scale(_h_br[ix], 1e4/ *_netap);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_m,_h_br[3];
    CounterPtr _netap;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1504943);

}
