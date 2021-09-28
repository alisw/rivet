#// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi- Xibar+ cross section
  class BESIII_2019_I1758883 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2019_I1758883);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nXi, "/TMP/nXi" );
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
      // find the Xis
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	if(abs(p1.pid())!=3312) continue;
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
	  if(abs(p2.pid())!=3312) continue;
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
	    _nXi->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nXi->val();
      double error = _nXi->err();
      sigma *= crossSection()/ sumOfWeights() /femtobarn;
      error *= crossSection()/ sumOfWeights() /femtobarn; 
      Scatter2D temphisto(refData(1, 1, 5));
      Scatter2DPtr  mult;
      book(mult,1, 1, 5);
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

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nXi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESIII_2019_I1758883);


}
