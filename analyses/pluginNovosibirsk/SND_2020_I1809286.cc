// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi-pi0
  class SND_2020_I1809286 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2020_I1809286);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_c_total, "/TMP/total");
      book(_c_omega, "/TMP/omega");
      book(_c_rho  , "/TMP/rho"  );
      book(_c_rhop , "/TMP/rhop" );
      if(inRange(sqrtS()/GeV,1.42,1.48)) {
	book(_h_x,4,1,1);
	book(_h_m,4,1,3);
      }
      else if(inRange(sqrtS()/GeV,1.65,1.68)) {
	book(_h_x,4,1,2);
	book(_h_m,4,1,4);
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
      Particle pip,pim;
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(p.pid()     == 211) pip=p;
	else if(p.pid()==-211) pim=p;
      }
      if(ntotal!=3) vetoEvent;
      if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==1)
	_c_total->fill();
      else
	vetoEvent;
      if(_h_x) {
	_h_x->fill(pip.momentum().p()/sqrtS());
	_h_x->fill(pim.momentum().p()/sqrtS());
	_h_m->fill((pip.momentum()+pim.momentum()).mass()/MeV);
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==223 or
					     Cuts::abspid==113    or Cuts::abspid==213 or
					     Cuts::abspid==100113 or Cuts::abspid==100213)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=1) continue;
	int idOther = 211;
	if(p.pid()==223 || p.pid()==113 || p.pid()==100113)
	  idOther = 111;
	else if(p.pid()==213 || p.pid()==100213)
	  idOther = -211;
	bool matched=true;
	for(auto const & val : nRes) {
	  if(val.first==idOther ) {
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
	if(matched) {
	  if(p.pid()==223)
	    _c_omega->fill();
	  else if(p.pid()==213 || p.pid()==113)
	    _c_rho->fill();
	  else
	    _c_rhop->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_x) {
	normalize(_h_x,1.,false);
	normalize(_h_m,1.,false);
      }
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=1;ix<4;++ix) {
	unsigned int ymax = ix!=3 ? 2 : 4;
	for(unsigned int iy=1;iy<ymax;++iy) {
	  double sigma(0.),error(0.);
	  if(ix<3) {
	    sigma = _c_total->val()*fact;
	    error = _c_total->err()*fact;
	  }
	  else if(ix==3) {
	    if(iy==1) {
	      sigma = _c_rho->val()*fact;
	      error = _c_rho->err()*fact;
	    }
	    else if(iy==2) {
	      sigma = _c_rhop->val()*fact;
	      error = _c_rhop->err()*fact;
	    }
	    else {
	      sigma = _c_omega->val()*fact;
	      error = _c_omega->err()*fact;
	    }
	  }
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr  mult;
	  book(mult, ix, 1, iy);
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
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_total,_c_omega,_c_rho,_c_rhop;
    Histo1DPtr _h_x,_h_m;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(SND_2020_I1809286);

}
