// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKI_1975_I100592 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MARKI_1975_I100592);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_nHadrons, "TMP/hadrons");
      book(_nEvent, "TMP/event");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _nEvent->fill();
      _nHadrons->fill(fs.particles().size());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
	double sigma,error;
	if(ix==1) {
	  sigma = _nEvent->val()*crossSection()/ sumOfWeights() /nanobarn;
	  error = _nEvent->err()*crossSection()/ sumOfWeights() /nanobarn;
	}
	else {
	  sigma = _nHadrons->val()/sumOfWeights();
	  error = _nHadrons->err()/sumOfWeights();
	}
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr mult;
	book(mult, ix, 1, 1);
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
    CounterPtr _nHadrons,_nEvent;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MARKI_1975_I100592);


}
