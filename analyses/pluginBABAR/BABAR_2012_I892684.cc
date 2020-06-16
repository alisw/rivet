// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2012_I892684 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2012_I892684);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_cKpKmpippim , "TMP/KpKmpippim");
      book(_cKstarKpi   , "TMP/KstarKpi");
      book(_cphipippim  , "TMP/phipippim");
      book(_cphif0_980  , "TMP/phif0_980");
      book(_cphif0_600  , "TMP/phif0_600");
      book(_cKpKmpi0pi0 , "TMP/KpKmpi0pi0");
      book(_cphif0pi0pi0, "TMP/phif0pi0pi0");
      book(_c2Kp2Km     , "TMP/2Kp2Km");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// K*0
	if(abs(p.pid())==313) {
	  map<long,int> nRes=nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // K* K+/- pi-/+
	  if(ncount !=2 ) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==321 || abs(val.first)==211) {
	      continue;
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched==false) continue;
	  if((nCount[321] == 1 && nCount[-321] ==0 &&
	      nCount[211] == 0 && nCount[-211] == 1) ||
	     (nCount[321] == 0 && nCount[-321] ==1 &&
	      nCount[211] == 1 && nCount[-211] == 0))
	    _cKstarKpi->fill();
	}
	else if(p.pid()==333) {
	  map<long,int> nRes=nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // phi pi+pi-
	  if(ncount==2) {
	    bool matched = true;
	    for(auto const & val : nRes) {
	      if(abs(val.first)==211) {
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
	      _cphipippim->fill();
	  }
	  for (const Particle& p2 : ufs.particles()) {
	    if(p2.pid()!=9010221&&p2.pid()!=9000221) continue;
	    if(p2.parents()[0].isSame(p)) continue;
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
	      if(p2.pid()==9010221) {
		_cphif0pi0pi0->fill();
		_cphif0_980  ->fill();
	      }
	      else {
		_cphif0_600  ->fill();
	      }
	    }
	  }
	}
      }
      if(ntotal==4) {
	if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1)
	  _cKpKmpippim->fill();
	else if( nCount[321]==1 && nCount[-321]==1 && nCount[111]==2)
	  _cKpKmpi0pi0->fill();
	else if( nCount[321]==2 && nCount[-321]==2)
	  _c2Kp2Km->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<9;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma = _cKpKmpippim->val();
	  error = _cKpKmpippim->err();
	}
	else if(ix==2) {
	  sigma = _cKstarKpi->val();
	  error = _cKstarKpi->err();
	}
     	else if(ix==3) {
	  sigma = _cphipippim->val();
	  error = _cphipippim->err();
	}
     	else if(ix==4) {
	  sigma = _cphif0_980->val();
	  error = _cphif0_980->err();
	}
     	else if(ix==5) {
	  sigma = _cphif0_600->val();
	  error = _cphif0_600->err();
	}
     	else if(ix==6) {
	  sigma = _cKpKmpi0pi0->val();
	  error = _cKpKmpi0pi0->err();
	}
     	else if(ix==7) {
	  sigma = _cphif0pi0pi0->val();
	  error = _cphif0pi0pi0->err();
	}
     	else if(ix==8) {
	  sigma =  _c2Kp2Km->val();
	  error =  _c2Kp2Km->err();
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
    CounterPtr _cKpKmpippim, _cKstarKpi, _cphipippim,
      _cphif0_980,_cphif0_600, _cKpKmpi0pi0, _cphif0pi0pi0, _c2Kp2Km;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2012_I892684);


}
