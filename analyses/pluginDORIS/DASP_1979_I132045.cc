// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DASP_1979_I132045 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DASP_1979_I132045);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // find the hists based on beam energies
      int ihist=-1;
      if (inRange(sqrtS()/GeV,3.6,3.67)) {
	ihist=0;
      }
      else if (inRange(sqrtS()/GeV,3.98,4.1)) {
	ihist=1;
      }
      else if (inRange(sqrtS()/GeV,4.1,4.24)) {
	ihist=2;
      }
      else if (inRange(sqrtS()/GeV,4.24,4.36)) {
	ihist=3;
      }
      else if (inRange(sqrtS()/GeV,4.36,4.46)) {
	ihist=4;
      }
      else if (inRange(sqrtS()/GeV,4.46,4.98)) {
	ihist=5;
      }
      else if (isCompatibleWithSqrtS(5.0)) {
	ihist=6;
      }
      else if (isCompatibleWithSqrtS(5.2)) {
	ihist=7;
      }
      else {
	MSG_ERROR("Beam energy not supported!");
      }
      // Book histograms
      book(_h_pi_p    ,  1,1,1+ihist);
      book(_h_K_p     ,  2+ihist,1,1);
      book(_h_proton_p, 10+ihist,1,1);
      book(_h_pi_x    , 18,1,1+ihist);
      book(_h_K_x     , 19+ihist,1,1);
      book(_h_proton_x, 27+ihist,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
	const int id = p.abspid();
	const double modp = p.p3().mod();
	const double xp = 2.*modp/sqrtS();
	const double beta = modp / p.E();
	if(id==211) {
	  _h_pi_p->fill(modp);
	  _h_pi_x->fill(xp  ,1./beta);
	}
	else if(id==321) {
	  _h_K_p->fill(modp);
	  _h_K_x->fill(xp  ,1./beta);
	}
	else if(id==2212) {
	  _h_proton_p->fill(modp);
	  _h_proton_x->fill(xp  ,1./beta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      scale(_h_pi_p     ,fact);
      scale(_h_K_p      ,fact);
      scale(_h_proton_p ,fact);
      fact *= sqr(sqrtS());
      scale(_h_pi_x     ,fact);
      scale(_h_K_x      ,fact);
      scale(_h_proton_x ,fact);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pi_p, _h_K_p, _h_proton_p;
    Histo1DPtr _h_pi_x, _h_K_x, _h_proton_x;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DASP_1979_I132045);


}
