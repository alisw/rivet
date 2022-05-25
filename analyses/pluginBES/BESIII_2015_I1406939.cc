// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2015_I1406939 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1406939);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nChi1, "TMP/chi1");
      book(_nChi2, "TMP/chi2");
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
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool matched = false;
      for (const Particle& p :  ufs.particles(Cuts::pid==20443 or Cuts::pid==445)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	for (const Particle& p2 :  ufs.particles(Cuts::pid==223)) {
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(p.pid()==20443)
	      _nChi1->fill();
	    else if(p.pid()==445)
	      _nChi2->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact =  crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=2;ix<4;++ix) {
	double sigma(0.),error(0.);
	if(ix==2) {
	  sigma = _nChi1->val()*fact;
	  error = _nChi1->err()*fact;
	}
	else if(ix==3) {
	  sigma = _nChi2->val()*fact;
	  error = _nChi2->err()*fact;
	}
	Scatter2D temphisto(refData(ix, 1, 6));
	Scatter2DPtr  mult;
	book(mult, ix, 1, 6);
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
    CounterPtr _nChi1,_nChi2;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2015_I1406939);


}
