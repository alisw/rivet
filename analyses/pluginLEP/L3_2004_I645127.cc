// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief mu+mu- and tau+tau- in gamma gamma
  class L3_2004_I645127 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_2004_I645127);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // get the mode and options
      _mode =0;
      if( getOption("PROCESS") == "EE" ) _mode = 0;
      else if( getOption("PROCESS") == "GG") _mode = 1; 
      
      // Initialise and register projections
      if(_mode==0) {
	const GammaGammaKinematics& diskin = declare(GammaGammaKinematics(), "Kinematics");
	declare(GammaGammaFinalState(diskin), "FS");
	declare(UnstableParticles(),"UFS");
      }
      else if(_mode==1) {
	declare(FinalState(), "FS");
      }
      // Book counters
      book(_c_sigma_mu1, "/TMP/sigma_mu_1");
      book(_c_sigma_mu2, "/TMP/sigma_mu_2");
      book(_c_sigma_tau, "/TMP/sigma_tau");

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
      // stuff for e+e- collisions
      double W2 = sqr(sqrtS());
      if(_mode==0) {
	const GammaGammaKinematics& kin = apply<GammaGammaKinematics>(event, "Kinematics");
	W2 = kin.W2();
	if(W2<9.*sqr(GeV) ) vetoEvent;
      }
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      bool fiducal = true;
      for(const Particle & p : fs.particles()) {
     	nCount[p.pid()] += 1;
     	++ntotal;
	if(abs(p.pid())==13) {
	  if(abs(cos(p.momentum().polarAngle()))>0.8) fiducal = false;
	}
      }
      if( nCount[-13]==1 && nCount[13]==1 && ntotal==2+nCount[22]) {
	if(W2<1600.*sqr(GeV)) { 
	  _c_sigma_mu1->fill();
	  if(fiducal) {
	    _c_sigma_mu2->fill();
	  }
	} 
      }
      if(_mode==1) return;
      bool foundTauPlus = false, foundTauMinus = true;
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
    	// find the taus
     	if(abs(p.pid())==15) {
	  if(p.pid()== 15) foundTauMinus=true;
	  if(p.pid()==-15) foundTauPlus =true;
     	  findChildren(p,nCount,ntotal);
	}
      }
      if(!foundTauPlus || !foundTauMinus)
	vetoEvent;
      bool matched = true;
      for(auto const & val : nCount) {
	if(val.first==22) {
	  continue;
	}
	else if(val.second!=0) {
	  matched = false;
	  break;
	}
      }
      if(matched)
	_c_sigma_tau->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // prefactor for the cross sections
      double fact  = crossSection()/picobarn/sumOfWeights();
      if(_mode==1) fact /= 1000.;
      scale(_c_sigma_mu1,fact);
      scale(_c_sigma_mu2,fact);
      scale(_c_sigma_tau,fact);
      unsigned int imin=0, imax = 3;
      if(_mode==1) {
	imin=3;
	imax=4;
      }
      for(unsigned int ihist=imin;ihist<imax;++ihist) {
	unsigned int id=0, iy=0;
	double sigma = 0., error = 0.;
	if(ihist==0) {
	  id=1;
	  iy=1;
	  sigma = _c_sigma_mu2->val();
	  error = _c_sigma_mu2->err();
	}
	else if(ihist==1) {
	  id=1;
	  iy=2;
	  sigma = _c_sigma_mu1->val();
	  error = _c_sigma_mu1->err();
	}
	else if(ihist==2) {
	  id=2;
	  iy=1;
	  sigma = _c_sigma_tau->val();
	  error = _c_sigma_tau->err();
	}
	else if(ihist==3) {
	  id=3;
	  iy=5;
	  sigma = _c_sigma_mu1->val();
	  error = _c_sigma_mu1->err();
	}
	Scatter2D temphisto(refData(id, 1, iy));
	Scatter2DPtr  mult;
	book(mult, id, 1, iy);
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
    CounterPtr _c_sigma_mu1,_c_sigma_mu2,_c_sigma_tau;
    unsigned int _mode;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(L3_2004_I645127);


}
