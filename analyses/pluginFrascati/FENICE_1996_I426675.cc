// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class FENICE_1996_I426675 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(FENICE_1996_I426675);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	vetoEvent;
      _c_hadrons->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr hadrons;
      book(hadrons, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
       	const double x  = temphisto.point(b).x();
       	pair<double,double> ex = temphisto.point(b).xErrs();
       	pair<double,double> ex2 = ex;
       	if(ex2.first ==0.) ex2. first=0.0001;
       	if(ex2.second==0.) ex2.second=0.0001;
      	if (inRange(sqr(sqrtS()/GeV), x-ex2.first, x+ex2.second)) {
       	  hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
      	}
      	else {
       	  hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_hadrons;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(FENICE_1996_I426675);


}
