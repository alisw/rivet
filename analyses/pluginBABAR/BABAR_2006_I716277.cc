// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I716277 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I716277);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numEtaGamma, "TMP/EtaGamma");
      book(_numEtaPrimeGamma, "TMP/EtaPrimeGamma");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child: p.children()) {
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
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the eta/eta prime
        if(p.pid()!=221 && p.pid()!=331) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta pi+pi-
	if(ncount!=1) continue;
	bool matched = true;
	for(auto const & val : nRes) {
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
	  if(p.pid()==221)
	    _numEtaGamma->fill();
	  else
	    _numEtaPrimeGamma->fill();
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _numEtaGamma->val();
	  error = _numEtaGamma->err();
	}
	else {
	  sigma = _numEtaPrimeGamma->val();
	  error = _numEtaPrimeGamma->err();
	}
	sigma *= crossSection()/ sumOfWeights() /femtobarn;
	error *= crossSection()/ sumOfWeights() /femtobarn; 
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
    CounterPtr _numEtaGamma,_numEtaPrimeGamma;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I716277);


}
