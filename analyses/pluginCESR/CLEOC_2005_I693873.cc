// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEOC_2005_I693873 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2005_I693873);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npipi, "TMP/npipi");
      book(_nKK, "TMP/nKK");
      book(_nppbar, "TMP/nppbar");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal!=2) vetoEvent;
      
      if(nCount[211]==1 && nCount[-211]==1)
	_npipi->fill();
      else if(nCount[321]==1 && nCount[-321]==1)
	_nKK->fill();
      else if(nCount[2212]==1 && nCount[-2212]==1)
	_nppbar->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=1;ix<4;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma =  _npipi->val();
	  error =  _npipi->err();
	}
	else if(ix==2) {
	  sigma = _nKK->val();
	  error = _nKK->err();
	}
     	else if(ix==3) {
	  sigma = _nppbar->val();
	  error = _nppbar->err();
	}
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn;
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr mult;
	book(mult, 1, 1, ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }

    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _npipi,_nKK,_nppbar;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CLEOC_2005_I693873);


}
