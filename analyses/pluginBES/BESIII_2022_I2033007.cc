// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> K+K-pi0
  class BESIII_2022_I2033007 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2033007);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Book histograms
      for(unsigned int ix=0;ix<4;++ix) {
	stringstream ss;
	ss << "TMP/n" << ix;
	book(_nMeson[ix], ss.str());
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==3 && nCount[321]==1 &&
	 nCount[-321]==1 && nCount[111]==1)
	_nMeson[0]->fill();
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==323 or Cuts::abspid==325 or Cuts::pid==333)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=1) continue;
	// phi pi0
	long id = p.pid()>0 ? -321 : 321;
	if (p.pid()==333) id = 111;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(val.first==id) {
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
	  if(p.pid()==333)
	    _nMeson[1]->fill();
	  else if (p.abspid()==325)
	    _nMeson[2]->fill();
	  else
	    _nMeson[3]->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix) {
	double sigma = _nMeson[ix]->val();
	double error = _nMeson[ix]->err();
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn; 
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
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2033007);

}
