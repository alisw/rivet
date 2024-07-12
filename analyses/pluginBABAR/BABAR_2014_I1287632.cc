// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  eta_c -> K+K- eta or pi0
  class BABAR_2014_I1287632 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2014_I1287632);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 441);
      declare(ufs, "UFS");
      DecayedParticles ETAC(ufs);
      ETAC.addStable(PID::PI0);
      ETAC.addStable(PID::ETA);
      declare(ETAC,"ETAC");
      // histos
      book(_h_KK[0],1,1,1);
      book(_h_KK[1],2,1,1);
      book(_h_Kpeta,1,1,2);
      book(_h_Kmeta,1,1,3);
      book(_h_Kppi ,2,1,2);
      book(_h_Kmpi ,2,1,3);
      book(_dalitz[0], "dalitz_1",50,0.5,7.0,50,0.5 ,7.0);
      book(_dalitz[1], "dalitz_2",50,0.2,6.5,50,0.2,6.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 321,1},{-321,1}, {111,1}};
      static const map<PdgId,unsigned int> & mode2 = { { 321,1},{-321,1}, {221,1}};
      DecayedParticles ETAC = apply<DecayedParticles>(event, "ETAC");
      // loop over particles
      for(unsigned int ix=0;ix<ETAC.decaying().size();++ix) {
	int imode=0;
	if     (ETAC.modeMatches(ix,3,mode1)) imode=0;
	else if(ETAC.modeMatches(ix,3,mode2)) imode=1;
	else continue;
	const Particle & Kp = ETAC.decayProducts()[ix].at( 321)[0];
	const Particle & Km = ETAC.decayProducts()[ix].at(-321)[0];
	if(imode==0 && ETAC.decaying()[ix].mass()>2.922 && ETAC.decaying()[ix].mass()<3.036) {
	  const Particle & pi0 = ETAC.decayProducts()[ix].at( 111)[0];
	  double mplus  = (Kp.momentum()+pi0.momentum()).mass2();
	  double mminus = (Km.momentum()+pi0.momentum()).mass2();
	  double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	  _h_KK[1]->fill(mKK);
	  _h_Kppi->fill(mplus);
	  _h_Kmpi->fill(mminus);
	  _dalitz[1]->fill(mplus,mminus);
	}
	else if (imode==1 && ETAC.decaying()[ix].mass()>2.910 && ETAC.decaying()[ix].mass()<3.03) {
	  const Particle & eta = ETAC.decayProducts()[ix].at( 221)[0];
	  double mplus  = (Kp.momentum()+eta.momentum()).mass2();
	  double mminus = (Km.momentum()+eta.momentum()).mass2();
	  double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	  _h_KK[0]->fill(mKK);
	  _h_Kpeta->fill(mplus);
	  _h_Kmeta->fill(mminus);
	  _dalitz[0]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK[0]);
      normalize(_h_KK[1]);
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_h_Kmeta);
      normalize(_h_Kpeta);
      normalize(_dalitz[0]);
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK[2],_h_Kmeta,_h_Kpeta,_h_Kmpi,_h_Kppi;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2014_I1287632);

}
