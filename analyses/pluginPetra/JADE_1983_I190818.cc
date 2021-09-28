// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Average multiplcity at a range of energies
  class JADE_1983_I190818 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1983_I190818);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if( !(fuzzyEquals(sqrtS()/GeV,12.0) ||
	    fuzzyEquals(sqrtS()/GeV,30.0) ||
	    fuzzyEquals(sqrtS()/GeV,35.0) )) {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book(_counter, "/TMP/MULT");
      book(_mult, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _counter->fill(cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_counter,1./sumOfWeights());

      double val = _counter->val();
      double err = _counter->err();
      
      Scatter2D tempScat(refData(1, 1, 1));
      
      for (size_t b = 0; b < tempScat.numPoints(); b++) {
        const double x  = tempScat.point(b).x();
        pair<double,double> ex = tempScat.point(b).xErrs();
        pair<double,double> ex2 = ex;
        if(ex2.first ==0.) ex2. first=0.0001;
        if(ex2.second==0.) ex2.second=0.0001;
        if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
          _mult->addPoint(x, val, ex, make_pair(err,err));
        }
        else {
          _mult->addPoint(x, 0., ex, make_pair(0.,.0));
        }
      }
    }
    //@}

  private:

    // Histogram
    CounterPtr _counter;
    Scatter2DPtr _mult;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1983_I190818);

}
