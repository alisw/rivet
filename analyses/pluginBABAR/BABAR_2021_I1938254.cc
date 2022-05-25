// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi-4pi0 and pi+pi-3pi0 eta
  class BABAR_2021_I1938254 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2021_I1938254);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Histograms
      for(unsigned int ix=0;ix<5;++ix) {
	book(_num[ix],"TMP/num_"+to_string(ix));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
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
      // stable final state
      if(nCount[211]==1 && nCount[-211]==1 && nCount[111]==4)
	_num[0]->fill();
      // intermediate states
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==221 || Cuts::pid==223)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta +X
	int idOther;
	if(p.pid()==221) {
	  idOther=223;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211 || val.first==-211 ) {
	      if(val.second !=1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.first==111) {
	      // 1 or 3 pi0
	      if(val.second !=3 && val.second != 1) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(nRes[111]==1)
	      _num[1]->fill();
	    else
	      _num[4]->fill();
	  }


	}
	// omega+X
	else {
	  idOther=221;
	  // omega+ 3pi0
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111) {
	      if(val.second !=3) {
		matched = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    _num[3]->fill();
	  }
	}
	for (const Particle& p2 : ufs.particles(Cuts::pid==idOther)) {
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=0) continue;
	  bool matched2 = true;
	  for(auto const & val : nResB) {
	    if(val.second!=0) {
	      matched2 = false;
	      break;
	    }
	  }
	  if(matched2) {
	    _num[2]->fill();
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<5;++ix) {
	double sigma = _num[ix]->val()*fact;
	double error = _num[ix]->err()*fact;
	Scatter2D temphisto(refData(1+ix, 1, 1));
	Scatter2DPtr  mult;
	book(mult, 1+ix, 1, 1);
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _num[5];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2021_I1938254);

}
