// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-
  class BELLE_2007_I749358 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I749358);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if(sqrtS()<0.8*GeV || sqrtS()>1.5*GeV)
	throw Error("Invalid CMS energy for BELLE_2007_I749358");
      // bin for the angle plots
      int ibin = (sqrtS()-0.8)/0.005 + 2;
      book(_h_cTheta,ibin,1,1);
      book(_cPi, "/TMP/nPi");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for(const Particle & p : part) {
	if(p.pid()==PID::PIPLUS) {
	  foundP=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
	else if(p.pid()==PID::PIMINUS)
	  foundM=true;
      }
      if(!foundP || !foundM) vetoEvent;
      if(cTheta<=0.6)    _cPi->fill();
      if(_h_cTheta ) _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_cTheta ) scale(_h_cTheta ,fact);
      double sigma = _cPi->val()*fact;
      double error = _cPi->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	  mult->addPoint(x, sigma, ex, make_pair(error,error));
	}
	else {
	  mult->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_cTheta;
    CounterPtr _cPi;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I749358);

}
