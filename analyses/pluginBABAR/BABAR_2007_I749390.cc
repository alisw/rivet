// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K+K-pi0
  class BABAR_2007_I749390 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I749390);


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
      book(_h_Kppi ,1,1,1);
      book(_h_Kmpi ,1,1,2);
      book(_h_KK   ,1,1,3);
      book(_dalitz, "dalitz",50,0.3,2.1,50,0.3,2.1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-321,1}, {111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( ! D0.modeMatches(ix,3,mode)  ) continue;
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & pi0 = D0.decayProducts()[ix].at(111);
	const Particles & Kp  = D0.decayProducts()[ix].at( sign*321);
	const Particles & Km  = D0.decayProducts()[ix].at(-sign*321);
	double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	_h_KK  ->fill(mKK);
	_h_Kppi->fill(mplus);
	_h_Kmpi->fill(mminus);
	_dalitz->fill(mminus,mplus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK  );
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK,_h_Kmpi,_h_Kppi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I749390);

}
