// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class GAMMAGAMMA_1981_I158474 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(GAMMAGAMMA_1981_I158474);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_n3pi, "TMP/n3pi");
      book(_n4pi, "TMP/n4pi");
      book(_n5pi, "TMP/n5pi");
      book(_n6pi, "TMP/n6pi");
      book(_n35pi, "TMP/n35pi");
      book(_n46pi, "TMP/n46pi");
      book(_nC2, "TMP/nC2");
      book(_nC4, "TMP/nC4");
      book(_nmu, "TMP/nmu");
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
	_nmu->fill();
      else {
	if(ntotal==3 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==1 ) {
	  _n3pi->fill();
	}
	if(ntotal==4 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==2 ) {
	  _n4pi->fill();
	}
	if(ntotal==5 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==1 ) {
	  _n5pi->fill();
	}
	if(ntotal==6 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==2 ) {
	  _n6pi->fill();
	}
	if(nCount[211] == 1 && nCount[-211]==1 && ntotal == 2+nCount[111]) {
	  _nC2->fill();
	}
	if(nCount[211] == 2 && nCount[-211]==2 && ntotal == 4+nCount[111]) {
	  _nC4->fill();
	}
	if((nCount[211]+nCount[-211]+nCount[111])==ntotal ) {
	  if(ntotal==3 || ntotal ==5)
	    _n35pi->fill();
	  else if(ntotal==4 || ntotal==6) 
	    _n46pi ->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=1;ix<7;++ix) {
	double sigma = 0.0, error = 0.0;
	if(ix==1) {
	  sigma = _n3pi->val()*fact;
	  error = _n3pi->err()*fact;
	}
	else if(ix==2) {
	  sigma = _n4pi->val()*fact;
	  error = _n4pi->err()*fact;
	}
	else if(ix==3) {
	  sigma = _n5pi->val()*fact;
	  error = _n5pi->err()*fact;
	}
	else if(ix==4) {
	  sigma = _n6pi->val()*fact;
	  error = _n6pi->err()*fact;
	}
	else if(ix==5) {
	  sigma = _n35pi->val()*fact;
	  error = _n35pi->err()*fact;
	}
	else if(ix==6) {
	  sigma = _n46pi->val()*fact;
	  error = _n46pi->err()*fact;
	} 
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
      for(unsigned int ix=1;ix<3;++ix) {
	Scatter1D R = (ix==1? *_nC2 : *_nC4)/ *_nmu;
	double              rval = R.point(0).x();
	pair<double,double> rerr = R.point(0).xErrs();
	double sig_h = (ix ==1 ? _nC2 : _nC4)->val()*fact;
	double err_h = (ix ==1 ? _nC2 : _nC4)->err()*fact;
	double sig_m = _nmu->val()*fact;
	double err_m = _nmu->err()*fact;
	Scatter2D temphisto(refData(2, 1, ix));
	std::ostringstream title;
	if(ix==1)
	  title << "sigma_2pi";
	else
	  title << "sigma_4pi";
	Scatter2DPtr hadrons;
	book(hadrons, title.str());
	Scatter2DPtr muons;
 book(muons, "sigma_muons");
	Scatter2DPtr mult;
	book(mult, 2,1,ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult   ->addPoint(x, rval, ex, rerr);
	    hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
	    if(ix==1) muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	  }
	  else {
	    mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	    if(ix==1) muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _n3pi,_n4pi,_n5pi,_n6pi,_n35pi,_n46pi,_nC2,_nC4,_nmu;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(GAMMAGAMMA_1981_I158474);


}
