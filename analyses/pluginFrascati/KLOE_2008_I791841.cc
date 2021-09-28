// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class KLOE_2008_I791841 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(KLOE_2008_I791841);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_n4pi, "TMP/4pi");
      book(_n2pigamma, "TMP/2pigamma");
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
      if(nCount[111]==2) {
	if( nCount[211] == 1 && nCount[-211] == 1 )
	  _n4pi->fill();
	else if( nCount[22] == 1)
	  _n2pigamma->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma = _n4pi->val();
	  error = _n4pi->err();
	}
	else if(ix==2) {
   	  sigma = _n2pigamma->val();
	  error = _n2pigamma->err();
        }
	sigma *= crossSection()/ sumOfWeights() /nanobarn;
	error *= crossSection()/ sumOfWeights() /nanobarn; 
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr mult;
	book(mult, ix, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
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
    CounterPtr _n4pi,_n2pigamma;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(KLOE_2008_I791841);


}
