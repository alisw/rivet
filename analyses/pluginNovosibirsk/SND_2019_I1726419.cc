// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief SND pi+pi-pi0eta cross section
  class SND_2019_I1726419 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2019_I1726419);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");


      book(_c_all  , "/TMP/all");
      book(_c_omega, "/TMP/omega");
      book(_c_phi  , "/TMP/phi");
      book(_c_rho  , "/TMP/rho");
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
      bool found = false, foundOmegaPhi = false;
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta pi+pi-
	if(ncount==3) {
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211 || val.first==111 ) {
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
	  if(matched) {
	    _c_all->fill();
	    found = true;
	  }
	}
	for (const Particle& p2 : ufs.particles()) {
	  if(p2.pid()!=223 && p2.pid()!=333) continue;
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
	    if(p2.pid()==223)
	      _c_omega->fill();
	    else if(p2.pid()==333)
	      _c_phi->fill();
	    foundOmegaPhi=true;
	  }
	}
      }
      if(found && !foundOmegaPhi)
	_c_rho->fill();
    }
    

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=1;ix<5;++ix) {
	double sigma(0.),error(0.);
	if(ix==1) {
	  sigma = _c_all->val()*fact;
	  error = _c_all->err()*fact;
	}
	else if(ix==2) {
	  sigma = _c_omega->val()*fact;
	  error = _c_omega->err()*fact;
	}
	else if(ix==3) {
	  sigma = _c_phi->val()*fact;
	  error = _c_phi->err()*fact;
	}
	else if(ix==4) {
	  sigma = _c_rho->val()*fact;
	  error = _c_rho->err()*fact;
	} 
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr  mult;
	book(mult, 1, 1, ix);
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
    CounterPtr _c_all, _c_omega, _c_phi,_c_rho;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2019_I1726419);


}
