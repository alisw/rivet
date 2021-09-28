// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2000_I503946 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SND_2000_I503946);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numOmegaPi, "TMP/OmegaPi");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]+=1;
	  ++ncount;
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
      // three particles (pi0 pi0 gamma)
      if(ntotal!=3) vetoEvent;
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the omega
	if(p.pid()==223) {
	  map<long,int> nRes;
	  int ncount(0);
	  findChildren(p,nRes,ncount);
	  // only omega to pi0 gamma mode
	  if(ncount!=2) continue;
	  if(nRes[111]!=1 || nRes[22]!=1) continue;
	  // omega pi0
	  if(nCount[111]-nRes[111]==1)
	    _numOmegaPi->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sigma = _numOmegaPi->val();
      double error = _numOmegaPi->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn; 
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	  mult->addPoint(x, sigma, ex, make_pair(error,error));
	}
	else {
	  mult->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _numOmegaPi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SND_2000_I503946);


}
