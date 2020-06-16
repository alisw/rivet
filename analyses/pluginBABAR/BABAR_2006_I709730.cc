// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I709730 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2006_I709730);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_num3pip3pim,      "TMP/num3pip3pim"     );
      book(_num2pip2pim2pi0,  "TMP/num2pip2pim2pi0" );
      book(_num2pip2pim2KpKm, "TMP/num2pip2pim2KpKm");
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

      if(ntotal!=6) vetoEvent;
      if(nCount[-211]==3 && nCount[211]==3)
	_num3pip3pim->fill();
      else if(nCount[-211]==2 && nCount[211]==2 && nCount[111]==2)
	_num2pip2pim2pi0->fill();
      else if(nCount[-211]==2 && nCount[211]==2 && nCount[321]==1 && nCount[-321]==1)
	_num2pip2pim2KpKm->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=1; ix<4; ++ix) {
        double sigma = 0., error = 0.;
        if(ix==1) {
          sigma = _num3pip3pim->val();
          error = _num3pip3pim->err();
        }
        else if(ix==2) {
          sigma = _num2pip2pim2pi0->val();
          error = _num2pip2pim2pi0->err();
        }
        else if(ix==3) {
          sigma = _num2pip2pim2KpKm->val();
          error = _num2pip2pim2KpKm->err();
        }
        sigma *= crossSection()/ sumOfWeights() /nanobarn;
        error *= crossSection()/ sumOfWeights() /nanobarn;
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

    // just count the number of events of the types we're looking for
    CounterPtr _num3pip3pim,_num2pip2pim2pi0,_num2pip2pim2KpKm;
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2006_I709730);


}
