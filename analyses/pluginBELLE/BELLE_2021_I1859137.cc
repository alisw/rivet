// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B(*) bar{B}(*) exclusive cross section
  class BELLE_2021_I1859137 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1859137);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_nBB  , "/TMP/nBB"  );
      book(_nBBS , "/TMP/nBBS" );
      book(_nBSBS, "/TMP/nBSBS");
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // extract botton hadrons
      Particles bHadrons=apply<FinalState>(event, "UFS").particles(Cuts::abspid==511 or Cuts::abspid==513 or
								   Cuts::abspid==521 or Cuts::abspid==523);
      for(unsigned int ix=0;ix<bHadrons.size();++ix) {
	long pix = bHadrons[ix].parents()[0].abspid();
	if(pix==511 || pix==413 || pix==521 || pix==523) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(bHadrons[ix],nRes,ncount);
	bool matched=false;
	for(unsigned int iy=ix+1;iy<bHadrons.size();++iy) {
	  long piy = bHadrons[ix].parents()[0].abspid();
	  if(piy==511 || piy==413 || piy==521 || piy==523) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(bHadrons[iy],nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(bHadrons[ix].abspid()==511 ||
	       bHadrons[ix].abspid()==521) {
	      if(bHadrons[iy].pid()==-bHadrons[ix].pid())
		_nBB->fill();
	      else
		_nBBS->fill();
	    }
	    else if(bHadrons[ix].abspid()==513 ||
		    bHadrons[ix].abspid()==523) {
	      if(bHadrons[iy].pid()==-bHadrons[ix].pid())
		_nBSBS->fill();
	      else
		_nBBS->fill();
	    }
	    break;
	  }
	}
	if(matched) break;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=1;iy<4;++iy) {
	double sigma,error;
	if(iy==1) {
	  sigma = _nBB->val();
	  error = _nBB->err();
	}
	else if(iy==2) {
	  sigma = _nBBS->val();
	  error = _nBBS->err();
	}
	else {
	  sigma = _nBSBS->val();
	  error = _nBSBS->err();
	}
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn; 
	Scatter2D temphisto(refData( 1, 1, iy));
    	Scatter2DPtr  mult;
        book(mult, 1, 1, iy);
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
    CounterPtr _nBB,_nBBS,_nBSBS;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1859137);

}
