// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > eta eta gamma cross section
  class SND_2021_I1942539 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2021_I1942539);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numEtaEtaGamma, "TMP/EtaEtaGamma");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
	 
      Particles etas = apply<FinalState>(event, "UFS").particles(Cuts::pid==221);
      // find the first eta
      for(unsigned int ix=0;ix<etas.size();++ix) {
	bool matched = false;
	if(etas[ix].children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(etas[ix],nRes,ncount);
	// find the second eta
	for(unsigned int iy=ix+1;iy<etas.size();++iy) {
	  if(etas[iy].children().empty()) continue;
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(etas[iy],nResB,ncountB);
	  matched = true;
	  for(auto const & val : nResB) {
	    if(val.first==22) {
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
	    _numEtaEtaGamma->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _numEtaEtaGamma->val();
      double error = _numEtaEtaGamma->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _numEtaEtaGamma;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(SND_2021_I1942539);

}
