// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-/K+ K-
  class ALEPH_2003_I626022 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2003_I626022);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // check CMS energy in range
      if(sqrtS()<2.*GeV || sqrtS()>6.*GeV)
	throw Error("Invalid CMS energy for ");
      // Final state
      declare(FinalState(),"FS");
      // histos
      book(_h_Pi,1,1,1);
      if(sqrtS()<4.)
	book(_h_K,2,1,1);
      book(_cPi, "/TMP/nPi_");
      book(_cK , "/TMP/nK_" );
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      if(part[0].pid()!=-part[1].pid()) vetoEvent;
      double cTheta(0.);
      bool foundPi(false),foundK(false);
      for(const Particle & p : part) {
	if(p.pid()==PID::PIPLUS) {
	  foundPi=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
	else if(p.pid()==PID::KPLUS) {
	  foundK=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
      }
      if(!foundPi && !foundK) vetoEvent;
      if(foundPi&&_h_Pi) _h_Pi->fill(cTheta);
      if(foundK &&_h_K ) _h_K ->fill(cTheta);
      if(foundPi)     _cPi->fill();
      else if(foundK)  _cK->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_Pi) scale(_h_Pi, fact);
      if(_h_K ) scale(_h_K , fact);
      for(unsigned int ih=3;ih<5;++ih) {
	CounterPtr count = ih==3 ? _cPi : _cK;
	double sigma = count->val()*fact;
	double error = count->err()*fact;
	// hist for axis
	Scatter2D temphisto(refData(ih, 1, 1));
	Scatter2DPtr cross;
	book(cross, ih, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	    cross->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    cross->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }
    
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Pi,_h_K;
    CounterPtr _cPi,_cK;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALEPH_2003_I626022);

}
