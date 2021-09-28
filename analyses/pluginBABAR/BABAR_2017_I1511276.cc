// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2017_I1511276 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2017_I1511276);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_nKKpi  , "TMP/KKpi");
      book(_nPhipi , "TMP/Phipi");
      book(_nKKeta , "TMP/KKeta");
      book(_nKKpipi, "TMP/KKpipi");

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
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
      // stable histos
      if( ntotal == 3 && nCount[130] == 1 &&
	  nCount[310] == 1 && nCount[111] == 1)
	_nKKpi->fill();
      else if( ntotal == 4 && nCount[130] == 1 &&
	       nCount[310] == 1 && nCount[111] == 2)
	_nKKpipi->fill();
      // unstable particles
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	if(p.pid()!=333 && p.pid()!=221) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	bool matched  = true;
	if(p.pid()==333 && ncount==1) {
	  for(auto const & val : nRes) {
	    if(val.first==111) {
	      if(val.second!=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched)
	    _nPhipi->fill();
	}
	else if(p.pid()==221 && ncount==2) {
	  for(auto const & val : nRes) {
	    if(val.first==130 || val.first==310) {
	      if(val.second!=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched)
	    _nKKeta->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<5;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma = _nKKpi->val();
	  error = _nKKpi->err();
	}
	else if (ix==2) {
	  sigma = _nPhipi->val();
	  error = _nPhipi->err();
	}
	else if (ix==3) {
	  sigma = _nKKeta->val();
	  error = _nKKeta->err();
	}
	else if (ix==4) {
	  sigma = _nKKpipi->val();
	  error = _nKKpipi->err();
	}
	  
	sigma *= crossSection()/ sumOfWeights() /nanobarn;
	error *= crossSection()/ sumOfWeights() /nanobarn;
	
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr  mult;
        book(mult, ix, 1, 1);
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
    CounterPtr _nKKpi,_nPhipi,_nKKeta,_nKKpipi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2017_I1511276);


}
