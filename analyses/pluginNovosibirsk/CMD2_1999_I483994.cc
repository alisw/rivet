// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD2_1999_I483994 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_1999_I483994);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_ncharged, "TMP/charged");
      book(_nneutral, "TMP/neutral");
      book(_nomega, "TMP/omega");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      if(ntotal==4) {
	if( nCount[211] == 2 && nCount[-211] == 2 )
	  _ncharged->fill();
	else if( nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 2)
	  _nneutral->fill();
      }
      

      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the omega
	if(p.pid()==223) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // omega pi+pi-
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111) {
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
	  if(matched)
	    _nomega->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<4;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma = _ncharged->val();
	  error = _ncharged->err();
	}
	else if(ix==2) {
   	  sigma = _nneutral->val();
	  error = _nneutral->err();
        }
	else if(ix==3) {
	  sigma = _nomega->val();
	  error = _nomega->err();
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
    CounterPtr _ncharged,_nneutral,_nomega;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMD2_1999_I483994);


}
