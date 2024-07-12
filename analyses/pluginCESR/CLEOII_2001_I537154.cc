// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D0 -> K- pi+ pi0
  class CLEOII_2001_I537154 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_2001_I537154);


    /// @name Analysis methods
    //@{

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
      book(_h_Kpi0,1,1,1);
      book(_h_Kpip,1,2,1);
      book(_h_pipi,1,3,1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double E0   = 22.1e-5  ;
      static const double Ex   = -6.89e-5 ;
      static const double Ey   = -27.1e-5 ; 
      static const double Ex2  = 10.4e-5  ;
      static const double Exy  = 38.2e-5  ;
      static const double Ey2  = 12.4e-5  ;
      static const double Ex3  = -3.00e-5;
      static const double Ex2y = -7.97e-5;
      static const double Exy2 = -12.8e-5;
      static const double Ey3  = -0.53e-5;

      static const map<PdgId,unsigned int> & mode   = { {-321,1},{ 211,1}, {111,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-211,1}, {111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.decaying()[ix].pid()>0 && D0.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (D0.decaying()[ix].pid()<0 && D0.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & pi0= D0.decayProducts()[ix].at(111);
	const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	double mKpip = (Km [0].momentum() + pip[0].momentum()).mass2();
	double mKpi0 = (Km [0].momentum() + pi0[0].momentum()).mass2();
	double mpipi = (pi0[0].momentum() + pip[0].momentum()).mass2();
	double eff = E0 + Ex*mKpip + Ey*mpipi +Ex2*sqr(mKpip)+Exy*mKpip*mpipi+Ey2*sqr(mpipi)
	  +Ex3*pow(mKpip,3)+Ex2y*sqr(mKpip)*mpipi +Exy2*mKpip*sqr(mpipi)+Ey3*pow(mpipi,3);
	_h_Kpi0->fill(mKpi0,eff);
	_h_Kpip->fill(mKpip,eff);
	_h_pipi->fill(mpipi,eff);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi0);
      normalize(_h_Kpip);
      normalize(_h_pipi);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Kpi0,_h_Kpip,_h_pipi;
    //@}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2001_I537154);

}
