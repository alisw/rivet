// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief  gamma gamma -> pi+pi-/K+ K-
  class BELLE_2005_I667712 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I667712);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if(sqrtS()<2.4*GeV || sqrtS()>4.1*GeV)
	throw Error("Invalid CMS energy for BELLE_2005_I667712");
      for(unsigned int ix=0;ix<7;++ix) {
	std::ostringstream title;
	title << "/TMP/nPi_" << ix;
	book(_cPi[ix], title.str());
      }
      for(unsigned int ix=0;ix<7;++ix) {
	std::ostringstream title;
	title << "/TMP/nK_" << ix;
	book(_cK[ix], title.str());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      if(part[0].pid()!=-part[1].pid()) vetoEvent;
      double cTheta(0.);
      bool foundPi(false),foundK(false);
      for(const Particle & p : part) {
	if(p.pid()==PID::PIPLUS) {
	  foundPi=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
	else if(p.pid()==PID::KPLUS) {
	  foundK=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
      }
      if(!foundPi && !foundK) vetoEvent;
      int ibin = cTheta/0.1;
      if(ibin>5) vetoEvent;
      if(foundPi) {
	_cPi[   0  ]->fill();
	_cPi[ibin+1]->fill();
      }
      else if(foundK) {
	_cK[   0  ]->fill();
	_cK[ibin+1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ip=0;ip<2;++ip) {
	CounterPtr denom = ip==0 ? _cPi[0] : _cK[0];
	if(denom->numEntries()==0) continue;
	for(unsigned int ih=0;ih<7;++ih) {
	  CounterPtr numer = ip==0 ? _cPi[ih] : _cK[ih];
	  double sigma = numer->val()*fact;
	  double error = numer->err()*fact;
	  // bin width for 2d dist
	  if(ih!=0) {
	    sigma /=0.1;
	    error /=0.1;
	  }
	  unsigned int ix=5+ip, iy=ih;
	  if(ih==0) {
	    ix=1;
	    iy = ip==0 ? 2 : 1;
	  }
	  // ratio
	  std::ostringstream title;
	  title << "/TMP/n_" << ip << "_" << ih;
	  Scatter1D sTemp(title.str());
	  Scatter1DPtr s1d = registerAO(sTemp);
	  // hist for axis
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr cross,ratio;
	  book(cross, ix, 1, iy);
	  if(ih!=0) {
	    book(ratio,ix-2,1,iy);
	    divide(numer,denom,s1d);
	  }
	  for (size_t b = 0; b < temphisto.numPoints(); b++) {
	    const double x  = temphisto.point(b).x();
	    pair<double,double> ex = temphisto.point(b).xErrs();
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	      cross->addPoint(x, sigma, ex, make_pair(error,error));
	      if(ih!=0) {
		ratio->addPoint(x,s1d->points()[0].x()/0.1,ex,make_pair(s1d->points()[0].xErrs().first /0.1,
									s1d->points()[0].xErrs().second/0.1));
	      }
	    }
	    else {
	      cross->addPoint(x, 0., ex, make_pair(0.,.0));
	      if(ih!=0)
		ratio->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	}
      }


      // }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPi[7],_cK[7];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I667712);

}
