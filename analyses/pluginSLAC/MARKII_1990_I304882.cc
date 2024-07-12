// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+ pi-
  class MARKII_1990_I304882 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1990_I304882);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if(sqrtS()<0.35*GeV || sqrtS()>1.6*GeV)
	throw Error("Invalid CMS energy for MARKII_1990_I304882");
      for(unsigned int ix=0;ix<7;++ix) {
	std::ostringstream title;
	title << "/TMP/nPi_" << ix;
	book(_cPi[ix], title.str());
      }
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
      int ibin = cTheta/0.1;
      if(ibin>5) vetoEvent;
      _cPi[   0  ]->fill();
      _cPi[ibin+1]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ih=0;ih<7;++ih) {
	unsigned int ix=2,iy=ih;
	if(ih==0) {
	  ix=1;
	  iy=1;
	}
	double sigma = _cPi[ih]->val()*fact;
	double error = _cPi[ih]->err()*fact;
	Scatter2D temphisto(refData(ix, 1, iy));
	Scatter2DPtr mult;
	book(mult, ix, 1, iy);
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
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPi[7];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1990_I304882);

}
