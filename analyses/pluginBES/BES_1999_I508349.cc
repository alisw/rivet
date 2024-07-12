// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief charm cross sections
  class BES_1999_I508349 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1999_I508349);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(),"UFS");
      book(_nD0   ,"/TMP/nD0"   );
      book(_nDp   ,"/TMP/nDp"   );
      book(_nDs   ,"/TMP/nDs"   );
      book(_nCharm,"/TMP/nCharm");
      if(isCompatibleWithSqrtS(4.03,1e-3)) {
	book(_h_D0,2,1,1);
	book(_h_Dp,2,1,2);
      }
      else if(isCompatibleWithSqrtS(4.14,1e-3)) {
	book(_h_D0,3,1,1);
	book(_h_Dp,3,1,2);
      }
      else
	MSG_ERROR("Beam energy not supported!");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double nD0(0),nDp(0),nDs(0);
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==411 or
										Cuts::abspid==421 or
										Cuts::abspid==431)) {
	if(p.abspid()==421) {
	  _h_D0->fill(p.momentum().p3().mod());
	  ++nD0;
	}
	else if(p.abspid()==411) {
	  _h_Dp->fill(p.momentum().p3().mod());
	  ++nDp;
	}
	else {
	  ++nDs;
	}
      }
      _nCharm->fill(0.5*(nD0+nDp+nDs));
      _nD0   ->fill(0.5*nD0);
      _nDp   ->fill(0.5*nDp);
      _nDs   ->fill(0.5*nDs);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_D0,0.5*crossSection()/sumOfWeights()/nanobarn);
      scale(_h_Dp,0.5*crossSection()/sumOfWeights()/nanobarn);
      for(unsigned int ix=1;ix<5;++ix) {
        double sigma = crossSection()/ sumOfWeights() /nanobarn;
        double error = crossSection()/ sumOfWeights() /nanobarn; 
      	if(ix==1) {
      	  sigma *= _nD0->val();
      	  error *= _nD0->err();
      	}
      	else if (ix==2){
      	  sigma *= _nDp->val();
      	  error *= _nDp->err();
      	}
      	else if (ix==3){
      	  sigma *= _nDs->val();
      	  error *= _nDs->err();
      	}
      	else if (ix==4){
      	  sigma *= _nCharm->val();
      	  error *= _nCharm->err();
      	}
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _nD0,_nDp,_nDs,_nCharm;
    Histo1DPtr _h_D0,_h_Dp;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BES_1999_I508349);

}
