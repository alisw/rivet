// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi-eta
  class BESIII_2022_I2039027 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2039027);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_num[0] , "TMP/etapipi");
      book(_num[1] , "TMP/etarho" );
    }
    
    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over eta mesons
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(!matched) continue;
	_num[0]->fill();
	for (const Particle& p2 : ufs.particles(Cuts::pid==113)) {
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=0) continue;
	  bool matched = true;
	  for(auto const & val : nResB) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) _num[1]->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<2;++ix) {
	double sigma = _num[ix]->val()*fact;
	double error = _num[ix]->err()*fact;
	Scatter2D temphisto(refData(ix+1, 1, 1));
	Scatter2DPtr  mult;
        book(mult, ix+1, 1, 1);
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

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _num[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2039027);

}
