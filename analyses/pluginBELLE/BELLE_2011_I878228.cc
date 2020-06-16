// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > Ds(*)+ Ds*()-
  class BELLE_2011_I878228 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2011_I878228);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(  _c_DpDm, "/TMP/sigma_DpDm"  );
      book( _c_DpDmS, "/TMP/sigma_DpDmS" );
      book(_c_DpSDmS, "/TMP/sigma_DpSDmS");
      book(   _c_All, "/TMP/sigma_All"   );
      book(    _c_mu, "/TMP/sigma_mu"    );
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
      	 ntotal==2+nCount[22]) {
	_c_mu->fill();
	return;
      }
      // unstable charm analysis
      
      Particles ds = apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==431 or Cuts::abspid==433);
      for(unsigned int ix=0;ix<ds.size();++ix) {
       	const Particle& p1 = ds[ix];
	int id1 = abs(p1.pid());
      	// check fs
      	bool fs = true;
      	for (const Particle & child : p1.children()) {
      	  if(child.pid()==p1.pid()) {
      	    fs = false;
      	    break;
      	  }
      	}
      	if(!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
      	bool matched=false;
       	int sign = p1.pid()/id1;
	// loop over the other fs particles
	for(unsigned int iy=ix+1;iy<ds.size();++iy) {
	  const Particle& p2 = ds[iy];
	  fs = true;
	  for (const Particle & child : p2.children()) {
	    if(child.pid()==p2.pid()) {
	      fs = false;
	      break;
	    }
	  }
	  if(!fs) continue;
	  if(p2.pid()/abs(p2.pid())==sign) continue;
	  int id2 = abs(p2.pid());
	  if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
	    continue;
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
	    if(id1==431 && id2==431) {
	      _c_DpDm->fill();
	    }
	    else if(id1==433 && id2==433) {
	      _c_DpSDmS->fill();
	    }
	    else if((id1==431 && id2==433) ||
		    (id1==433 && id2==431)) {
	      _c_DpDmS->fill();
	    }
	    _c_All->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/nanobarn;
      for(unsigned int iy=1;iy<4;++iy) {
	double sigma = 0.0, error = 0.0;
	if(iy==1) {
	  sigma = _c_DpDm->val()*fact;
	  error = _c_DpDm->err()*fact;
	}
	else if(iy==2) {
	  sigma = _c_DpDmS->val()*fact;
	  error = _c_DpDmS->err()*fact;
	}
	else if(iy==3) {
	  sigma = _c_DpSDmS->val()*fact;
	  error = _c_DpSDmS->err()*fact;
	}
	Scatter2D temphisto(refData(1, 1, iy));
        Scatter2DPtr     mult;
        book(mult, 1, 1, iy);
        for (size_t b = 0; b < temphisto.numPoints(); b++) {
          const double x  = temphisto.point(b).x();
          pair<double,double> ex = temphisto.point(b).xErrs();
          pair<double,double> ex2 = ex;
          if(ex2.first ==0.) ex2. first=0.0001;
          if(ex2.second==0.) ex2.second=0.0001;
          if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
            mult   ->addPoint(x, sigma, ex, make_pair(error,error));
          }
          else {
            mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
          }
        }
      }
      // R
      if(_c_mu->val()==0. || _c_All->val()==0.) return;
      Scatter1D R = *_c_All/ *_c_mu;
      double              rval = R.point(0).x();
      pair<double,double> rerr = R.point(0).xErrs();
      Scatter2D temphisto(refData(2, 1, 1));
      Scatter2DPtr mult;
      book(mult, 2,1,1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, rval, ex, rerr);
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_DpDm,_c_DpDmS,_c_DpSDmS,_c_All,_c_mu;
    //@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2011_I878228);

}
