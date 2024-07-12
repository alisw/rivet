// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > 2pi+ 2pi- 3pi0 and  2pi+ 2pi- 2pi0 eta
  class BABAR_2021_I1844422 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2021_I1844422);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_c_2pip2pim3pi0   ,"TMP/2pip2pim3pi0"   );
      book(_c_2pip2pimeta    ,"TMP/2pip2pimeta"    );
      book(_c_omegapi0eta    ,"TMP/omegapi0eta"    );
      book(_c_pippim2pi0omega,"TMP/pippim2pi0omega");
      book(_c_pippim2pi0eta  ,"TMP/pippim2pi0eta"  );
      book(_c_2pip2pim2pi0eta,"TMP/2pip2pim2pi0eta");
    }
    
    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      if(ntotal==7 && nCount[211]==2 && nCount[-211]==2 && nCount[111]==3)
	_c_2pip2pim3pi0->fill();
      // intermediate omega and eta mesons
      Particles unstable = apply<FinalState>(event, "UFS").particles(Cuts::pid==223 or Cuts::pid==221);
      for (const Particle& p : unstable) {
     	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta
	if(p.pid()==221) {
	  if(ncount==4) {
	    // 2pi+2pi- eta
	    if(nRes[211]==2 && nRes[-211]==2 ) {
	      _c_2pip2pimeta->fill();
	    }
	    // pi+pi- 2pi0 eta
	    else if(nRes[211]==1 && nRes[-211]==1 && nRes[111]==2 ) {
	      _c_pippim2pi0eta->fill();
	    }
	  }
	  else if(ncount==6) {
	    // 2pi+ 2pi- 2pi0 eta
	    if(nRes[211]==2 && nRes[-211]==2 && nRes[111]==2 ) {
	      _c_2pip2pim2pi0eta->fill();
	    }
	  }
	}
	// omega
	else if(p.pid()==223) {
	  // pi+pi- 2pi0 omega
	  if(ncount==4 && nRes[211]==1 && nRes[-211]==1 && nRes[111]==2 ) {
	    _c_pippim2pi0omega->fill();
	  }
	}
	// mode with both eta and omega
	// pi0 omega eta
	for (const Particle& p2 : unstable ) {
	  if(p2.pid()==p.pid()) continue;
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=1) continue;
	  bool matched = true;
	  for(auto const & val : nResB) {
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
	  if(matched) {
	    _c_omegapi0eta->fill();
	    break;
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<7;++ix) {
	double sigma = 0., error = 0.;
	if(ix==1) {
	  sigma = _c_2pip2pim3pi0->val();
	  error = _c_2pip2pim3pi0->err();
	}
	else if(ix==2) {
	  sigma = _c_2pip2pimeta->val();
	  error = _c_2pip2pimeta->err();
	}
	else if(ix==3) {
	  sigma = _c_omegapi0eta->val();
	  error = _c_omegapi0eta->err();
	}
	else if(ix==4) {
	  sigma = _c_pippim2pi0omega->val();
	  error = _c_pippim2pi0omega->err();
	}
	else if(ix==5) {
	  sigma = _c_pippim2pi0eta->val();
	  error = _c_pippim2pi0eta->err();
	}
	else if(ix==6) {
	  sigma = _c_2pip2pim2pi0eta->val();
	  error = _c_2pip2pim2pi0eta->err();
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_2pip2pim3pi0,  _c_2pip2pimeta, _c_omegapi0eta,
      _c_pippim2pi0omega, _c_pippim2pi0eta, _c_2pip2pim2pi0eta;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2021_I1844422);

}
