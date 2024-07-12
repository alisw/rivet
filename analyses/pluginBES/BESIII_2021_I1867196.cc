// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> D_s* D_sJ
  class BESIII_2021_I1867196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1867196);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_numD[ix],"TMP/num_"+to_string(ix));
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over D_s*
      for(const Particle & Dstar : ufs.particles(Cuts::abspid==433)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(Dstar,nRes,ncount);
	bool matched=false;
	for(const Particle & p : ufs.particles(Cuts::abspid==10431 or
					       Cuts::abspid==10433 or
					       Cuts::abspid==20433)) {
	  // check particle and antiparticle
	  if(Dstar.pid()*p.pid()>0) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p,nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(p.abspid()==10431) _numD[0]->fill();
	    else if(p.abspid()==20433) _numD[1]->fill();
	    else if(p.abspid()==10433) _numD[2]->fill();
	    break;
	  }
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	double sigma = _numD[ix]->val()*fact;
	double error = _numD[ix]->err()*fact;
	Scatter2D temphisto(refData(1+ix, 1, 1));
	Scatter2DPtr  mult;
	book(mult, 1+ix, 1, 1);
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
    CounterPtr _numD[3];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1867196);

}
