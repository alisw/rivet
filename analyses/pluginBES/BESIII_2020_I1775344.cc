// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > K+K- pi0pi0
  class BESIII_2020_I1775344 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1775344);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<6;++ix) {
	std::ostringstream title;
	title << "TMP/c_" << ix+1;
	book(_c[ix],title.str());
      }
      if(isCompatibleWithSqrtS(2.125,1e-3)) {
	book(_h_KK   ,7,1,1);
	book(_h_pipi ,7,1,2);
	book(_h_Kpi  ,7,1,3);
	book(_h_KKpi ,7,1,4);
	book(_h_Kpipi,7,1,5);
      }
      else if(isCompatibleWithSqrtS(2.396,1e-3)) {
	book(_h_KK   ,8,1,1);
	book(_h_pipi ,8,1,2);
	book(_h_Kpi  ,8,1,3);
	book(_h_KKpi ,8,1,4);
	book(_h_Kpipi,8,1,5);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      Particles Kp,pi0;
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(p.abspid()==321)
	  Kp.push_back(p);
	else if(p.pid()==111)
	  pi0.push_back(p);
      }
      // intermediates
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==100321 or
					     Cuts::abspid==10323  or
					     Cuts::abspid==20323  or
					     Cuts::pid   ==333    or
					     Cuts::abspid==323   )) {
	if(p.children().empty()) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// X-/+ with K+/-
	if((p.abspid()==100321 || p.abspid()== 10323 || p.abspid()==20323) && ncount==1) {
	  bool matched = true;
	  int Kid = -p.pid()/p.abspid()*321;
	  for(auto const & val : nRes) {
	    if(val.first==Kid) {
	      if(val.second!=1) {
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
	    if(p.abspid()==100321)
	      _c[2]->fill();
	    else if(p.abspid()== 20323)
	      _c[3]->fill();
	    else if(p.abspid()==10323)
	      _c[4]->fill();
	  }
	}
	// phi + 2pi0
	else if(p.pid()==333 && ncount==2) {
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(val.first==111) {
	      if(val.second!=2) {
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
	    _c[1]->fill();
	}
	// K*K*
	else if(p.abspid()==323) {
	  for (const Particle& p2 : ufs.particles(Cuts::pid==-p.pid())) {
	    map<long,int> nResB = nRes;
	    int ncountB = ncount;
	    findChildren(p2,nResB,ncountB);
	    if(ncountB!=0) continue;
	    bool matched = true;
	    for(auto const & val : nResB) {
	      if(val.second!=0) {
	    	matched = false;
	     	break;
	      }
	    }
	    if(matched)
	      _c[5]->fill();
	  }
	}
      }
      // final-state
      if(ntotal==4 && nCount[321]==1 && nCount[-321]==1 && nCount[111]==2) {
	_c[0]->fill();
	if(_h_KK) {
	  FourMomentum pKK = Kp[0].momentum()+Kp[1].momentum();
	  _h_KK->fill(pKK.mass());
	  FourMomentum pPi = pi0[0].momentum()+pi0[1].momentum();
	  _h_pipi->fill(pPi.mass());
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h_KKpi ->fill((pKK+pi0[ix].momentum()).mass());
	    _h_Kpipi->fill((pPi+ Kp[ix].momentum()).mass());
	    for(unsigned int iy=0;iy<2;++iy)
	      _h_Kpi->fill((Kp[ix].momentum()+pi0[iy].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix) {
	double sigma = _c[ix]->val()*crossSection()/ sumOfWeights() /nanobarn;;
	double error = _c[ix]->err()*crossSection()/ sumOfWeights() /nanobarn;;
	Scatter2D temphisto(refData(ix+1, 1, 1));
	Scatter2DPtr  mult;
        book(mult, ix+1, 1, 1);
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
      if(_h_KK) {
	normalize(_h_KK   );
	normalize(_h_pipi );
	normalize(_h_Kpi  );
	normalize(_h_KKpi );
	normalize(_h_Kpipi);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c[6];
    Histo1DPtr _h_KK, _h_pipi, _h_Kpi, _h_KKpi, _h_Kpipi;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1775344);

}
