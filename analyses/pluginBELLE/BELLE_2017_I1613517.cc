// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+/- D*-/+ and D*+ D*- cross sections
  class BELLE_2017_I1613517 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1613517);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_DpDmS, "/TMP/sigma_DpDmS");
      book(_c_DpSDmS, "/TMP/sigma_DpSDmS");
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
	 ntotal==2+nCount[22])
	vetoEvent;
      // unstable charm analysis
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
       	const Particle& p1 = ufs.particles()[ix];
       	int id1 = abs(p1.pid());
       	if(id1 != 411 && id1 != 413) continue;
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
      	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
      	  const Particle& p2 = ufs.particles()[iy];
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
       	  if(id2 != 411 && id2 != 413) continue;
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
	    if(id1==413 && id2==413) {
	      _c_DpSDmS->fill();
	    }
	    else if((id1==411 && id2==413) ||
		    (id1==413 && id2==411)) {
	      _c_DpDmS->fill();
	    }
	    break;
	  }
      	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/nanobarn;
      for(unsigned int iy=1;iy<3;++iy) {
	double sigma = 0.0, error = 0.0;
	if(iy==1) {
	  sigma = _c_DpDmS->val()*fact;
	  error = _c_DpDmS->err()*fact;
	}
	else if(iy==2) {
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
    }
    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_DpDmS, _c_DpSDmS;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BELLE_2017_I1613517);


}
