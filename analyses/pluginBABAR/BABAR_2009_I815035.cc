// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D D* and D* D* cross section
  class BABAR_2009_I815035 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I815035);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_D0_Dstar   ,"/TMP/c_D0_Dstar"   );
      book(_c_Dplus_Dstar,"/TMP/c_Dplus_Dstar");
      book(_c_D_Dstar    ,"/TMP/c_D_Dstar"    );
      book(_c_Dstar_Dstar,"/TMP/c_Dstar_Dstar");
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
      // unstable charm analysis
      Particles ds = apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==411 or Cuts::abspid==413 or
								      Cuts::abspid==421 or Cuts::abspid==423);
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
       	    if((id1==421 && id2==423) || (id1==423 && id2==421)) {
	      _c_D0_Dstar->fill();
	      _c_D_Dstar ->fill();
	    }
	    else if((id1==411 && id2==413) || (id1==413 && id2==411)) {
	      _c_Dplus_Dstar->fill();
	      _c_D_Dstar ->fill();
	    }
	    else if((id1==413 && id2==413) || (id1==423 && id2==423)) {
	      _c_Dstar_Dstar->fill();
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
      for(unsigned int ih=1;ih<5;++ih) {
	double sigma = 0.0, error = 0.0;
	unsigned int ix=1,iy=ih;
	if(ih==1) {
	  sigma = _c_D0_Dstar->val()*fact;
	  error = _c_D0_Dstar->err()*fact;
	}
	else if(ih==2) {
	  sigma = _c_Dplus_Dstar->val()*fact;
	  error = _c_Dplus_Dstar->err()*fact;
	}
	else if(ih==3) {
	  sigma = _c_D_Dstar->val()*fact;
	  error = _c_D_Dstar->err()*fact;
	}
	else if(ih==4) {
	  sigma = _c_Dstar_Dstar->val()*fact;
	  error = _c_Dstar_Dstar->err()*fact;
	  ix=2;
	  iy=1;
	}
	Scatter2D temphisto(refData(ix, 1, iy));
        Scatter2DPtr     mult;
        book(mult, ix, 1, iy);
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_D0_Dstar,_c_Dplus_Dstar,_c_D_Dstar, _c_Dstar_Dstar;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I815035);

}
