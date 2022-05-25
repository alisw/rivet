// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > sigma+- sigmabar -+
  class BESIII_2020_I1814783 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1814783);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_n_plus ,"/TMP/NPLUS" );
      book(_n_minus,"/TMP/NMINUS");
      if(isCompatibleWithSqrtS(2.396, 1E-2)) {
	book(_h_cTheta_A,3,1,1);
	book(_h_cTheta_B,3,1,2);
      }
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
      // find the Sigmas
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
      	const Particle& p1 = ufs.particles()[ix];
      	if(abs(p1.pid())!=3112&&abs(p1.pid())!=3222) continue;
      	bool matched = false;
       	// check fs
       	bool fs = true;
      	for(const Particle & child : p1.children()) {
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
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
       	  const Particle& p2 = ufs.particles()[iy];
       	  if(p2.pid() != -p1.pid()) continue;
       	  // check fs
      	  bool fs = true;
      	  for(const Particle & child : p2.children()) {
      	    if(child.pid()==p2.pid()) {
      	      fs = false;
      	      break;
      	    }
      	  }
      	  if(!fs) continue;
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
	    if(abs(p1.pid())==3222) {
	      _n_plus->fill();
	      if(_h_cTheta_A) {
		double cTheta = p1.pid()>0 ?
		  cos(p1.momentum().polarAngle()) :
		  cos(p2.momentum().polarAngle());
		_h_cTheta_A->fill(cTheta);
		_h_cTheta_B->fill(cTheta);
	      }
	    }
	    else if(abs(p1.pid())==3112)
	      _n_minus->fill();
	    break;
	  }
      	}
	if(matched) break;
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_cTheta_A) {
      	normalize(_h_cTheta_A);
      	normalize(_h_cTheta_B);
      }
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int iy=1;iy<3;++iy) {
	double sigma,error;
	if(iy==1) {
	  sigma = _n_plus->val()*fact;
	  error = _n_plus->err()*fact;
	}
	else {
	  sigma = _n_minus->val()*fact;
	  error = _n_minus->err()*fact;
	}
	for(unsigned int ix=1;ix<3;++ix) {
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
    CounterPtr _n_plus,_n_minus;
    Histo1DPtr _h_cTheta_A,_h_cTheta_B;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1814783);

}
