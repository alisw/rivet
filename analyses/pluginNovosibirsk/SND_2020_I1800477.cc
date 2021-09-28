// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > eta pi0 gamma
  class SND_2020_I1800477 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SND_2020_I1800477);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_c_omega, "/TMP/omega");
      book(_c_total, "/TMP/total");
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta
	if(p.pid()!=221) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// check eta pi0 gamma FS
	if(ncount!=2) continue;
	bool matched=true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==22 || val.first==111 ) {
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
	_c_total->fill();
	// check if from omega
	for (const Particle& p2 : ufs.particles()) {
	  if(p2.pid()!=223) continue;
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
	  if(matched2)
	    _c_omega->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=1;ix<4; ix+=2) {
       	double sigma(0.),error(0.);
       	if(ix==1) {
       	  sigma = _c_omega->val()*fact;
       	  error = _c_omega->err()*fact;
       	}
       	else if(ix==3) {
       	  sigma = _c_total->val()*fact;
       	  error = _c_total->err()*fact;
       	} 
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_total,_c_omega;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(SND_2020_I1800477);

}
