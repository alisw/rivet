// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K0 pi0pi0
  class CLEOC_2011_I913909 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2011_I913909);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      book(_h_Kpi ,1,1,1);
      book(_h_pipi,1,1,2);
      book(_dalitz, "dalitz",50,0.,1.9,50,0.3,3.1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {111,2}, { 310,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( ! D0.modeMatches(ix,3,mode  ) ) continue;
	const Particles & K0 = D0.decayProducts()[ix].at(310);
	const Particles & pi0= D0.decayProducts()[ix].at(111);
	double mpipi = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	_h_pipi->fill(mpipi);
	for(unsigned int ix=0;ix<2;++ix) {
	  double mKpi  = (K0[0].momentum()+pi0[ix].momentum()).mass2();
	  _h_Kpi->fill(mKpi);
	  _dalitz->fill(mpipi,mKpi);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi);
      normalize(_h_Kpi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi, _h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2011_I913909);

}
