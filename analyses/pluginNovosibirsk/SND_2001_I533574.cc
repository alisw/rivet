// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2001_I533574 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2001_I533574);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKpKm, "TMP/KpKm");
      book(_nK0K0, "TMP/K0K0");
      book(_n3pi, "TMP/3pi");
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
      if(ntotal==2) {
	if(nCount[321]==1 && nCount[-321]==1)
	  _nKpKm->fill();
	else if(nCount[130]==1 && nCount[310]==1)
	  _nK0K0->fill();
      }
      else if(ntotal==3 && nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 1)
	_n3pi->fill();

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=1;iy<3;++iy) {
	for(unsigned int ix=1;ix<5;++ix) {
	  double sigma = 0., error = 0.;
	  if(ix==1) {
	    sigma = _nKpKm->val();
	    error = _nKpKm->err();
	  }
	  else if(ix==2) {
	    sigma = _nK0K0->val();
	    error = _nK0K0->err();
	  }
	  else if(ix==3) {
	    sigma = _nK0K0->val();
	    error = _nK0K0->err();
	  }
	  else if(ix==4) {
	    sigma = _n3pi->val();
	    error = _n3pi->err();
	  }
	  sigma *= crossSection()/ sumOfWeights() /nanobarn;
	  error *= crossSection()/ sumOfWeights() /nanobarn; 
	  Scatter2D temphisto(refData(iy, 1, ix));
	  Scatter2DPtr mult;
	  book(mult, iy, 1, ix);
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
    }
    
    //@}


    /// @name Histograms
    //@{
    CounterPtr _nKpKm,_nK0K0,_n3pi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2001_I533574);


}
