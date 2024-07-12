// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi dalitz decays
  class BABAR_2017_I1512302 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1512302);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PI0);
      PSI.addStable(PID::K0S);
      declare(PSI,"PSI");
      // histos
      book(_h_pippim,1,1,1);
      book(_h_pippi0,1,1,2);
      book(_dalitz_3pi, "dalitz_3pi",50,0.,9.,50,0.0,9.);
      book(_h_KpKm ,2,1,1);
      book(_h_Kppi0,2,1,2);
      book(_dalitz_KpKmpi, "dalitz_KpKmpi",50,0.,7.,50,0.0,7.);
      book(_h_K0Kp ,3,1,1);
      book(_h_K0pip,3,1,2);
      book(_h_Kppip,3,1,3);
      book(_dalitz_K0Kppim, "dalitz_K0Kppim",50,0.,8.,50,0.,8.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 211,1},{-211,1}, {111,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1},{-321,1}, {111,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 321,1},{-211,1}, {310,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-321,1},{ 211,1}, {310,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if (PSI.modeMatches(ix,3,mode1)) {
	  const Particle & pip = PSI.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = PSI.decayProducts()[ix].at(-211)[0];
	  const Particle & pi0 = PSI.decayProducts()[ix].at( 111)[0];
	  double mminus = (pim.momentum()+pi0.momentum()).mass2();
	  double mplus  = (pip.momentum()+pi0.momentum()).mass2();
	  double mneut  = (pip.momentum()+pim.momentum()).mass2();
	  _h_pippim->fill(mneut );
	  _h_pippi0->fill(mplus );
	  _h_pippi0->fill(mminus);
	  _dalitz_3pi->fill(mplus,mminus);
	}
	else if (PSI.modeMatches(ix,3,mode2)) {
	  const Particle & Kp  = PSI.decayProducts()[ix].at( 321)[0];
	  const Particle & Km  = PSI.decayProducts()[ix].at(-321)[0];
	  const Particle & pi0 = PSI.decayProducts()[ix].at( 111)[0];
	  double mminus = (Km.momentum()+pi0.momentum()).mass2();
	  double mplus  = (Kp.momentum()+pi0.momentum()).mass2();
	  double mneut  = (Kp.momentum()+Km.momentum()).mass2();
	  _h_KpKm->fill(mneut );
	  _h_Kppi0->fill(mplus );
	  _h_Kppi0->fill(mminus);
	  _dalitz_KpKmpi->fill(mplus,mminus);
	}
	else {
	  int sign =1;
	  if     (PSI.modeMatches(ix,3,mode3  )) sign= 1;
	  else if(PSI.modeMatches(ix,3,mode3CC)) sign=-1;
	  else continue;
	  const Particle & Kp  = PSI.decayProducts()[ix].at( sign*321)[0];
	  const Particle & pim = PSI.decayProducts()[ix].at(-sign*211)[0];
	  const Particle & K0  = PSI.decayProducts()[ix].at(      310)[0];
	  double mplus  = (Kp.momentum()  +  K0.momentum()).mass2();
	  double mminus = (K0.momentum()  + pim.momentum()).mass2();
	  double mKK    = (Kp.momentum()  +  K0.momentum()).mass2();
	  _h_K0Kp ->fill(mKK);
	  _h_Kppip->fill(mplus);
	  _h_K0pip->fill(mminus);
	  _dalitz_K0Kppim->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pippim,1.,false);
      normalize(_h_pippi0,1.,false);
      normalize(_dalitz_3pi);
      normalize(_h_KpKm,1.,false);
      normalize(_h_Kppi0,1.,false);
      normalize(_dalitz_KpKmpi);
      normalize(_h_Kppip);
      normalize(_h_K0pip);
      normalize(_h_K0Kp );
      normalize(_dalitz_K0Kppim);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim,_h_pippi0;
    Histo2DPtr _dalitz_3pi;
    Histo1DPtr _h_KpKm,_h_Kppi0;
    Histo2DPtr _dalitz_KpKmpi;
    Histo1DPtr _h_Kppip,_h_K0pip,_h_K0Kp;
    Histo2DPtr _dalitz_K0Kppim;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1512302);

}
