// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi0 eta
  class CLEO_2004_I649917 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2004_I649917);


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
      book(_h_Kpi  ,1,1,1);
      book(_h_etapi,1,1,2);
      book(_h_Keta ,1,1,3);
      book(_dalitz, "dalitz",50,0.4,1.9,50,0.3,1.8);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { { 310,1}, { 221,1},{ 111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if ( ! D0.modeMatches(ix,3,mode)) continue;
	const Particles & pi0= D0.decayProducts()[ix].at(111);
	const Particles & eta= D0.decayProducts()[ix].at(221);
	const Particles & K0 = D0.decayProducts()[ix].at(310);
	double metapi = (eta[0].momentum()+pi0[0].momentum()).mass2();
	double metaK  = (eta[0].momentum()+ K0[0].momentum()).mass2();
	double mKpi   = ( K0[0].momentum()+pi0[0].momentum()).mass2();
	_dalitz ->fill(metapi,mKpi);
	_h_etapi->fill(metapi);
	_h_Keta ->fill(metaK );
	_h_Kpi  ->fill(mKpi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi  );
      normalize(_h_etapi);
      normalize(_h_Keta );
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi,_h_etapi,_h_Keta;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2004_I649917);

}
