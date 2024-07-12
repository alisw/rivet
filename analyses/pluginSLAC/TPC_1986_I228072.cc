// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-/K+ K-
  class TPC_1986_I228072 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TPC_1986_I228072);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if(sqrtS()<0.5*GeV || sqrtS()>3.5*GeV)
	throw Error("Invalid CMS energy for TPC_1986_I228072");
      book(_cPi[0],"/TMP/nPi_06");
      book(_cPi[1],"/TMP/nPi_03");
      book(_cK    ,"/TMP/nK"    );
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
      if(cTheta>0.6) vetoEvent;
      if(foundPi) {
	_cPi[0]->fill();
	if(cTheta<0.3)
	  _cPi[1]->fill();
      }
      else if(foundK)
	_cK->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ih=1;ih<6;++ih) {
	if(ih==2||ih==3) continue;
	double sigma(0.),error(0.);
	CounterPtr count;
	if(ih==1) {
	  count = _cPi[0];
	}
	else if(ih==4) {
	  count = _cPi[1];
	}
	else if(ih==5) {
	  count = _cK;
	}
	else
	  assert(false);
	if(count->numEntries()==0) continue;
	sigma = count->val()*fact;
	error = count->err()*fact;
	Scatter2DPtr cross;
	Scatter2D temphisto(refData(ih, 1, 1));
	book(cross, ih, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	    cross->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    cross->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPi[2],_cK;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(TPC_1986_I228072);

}
