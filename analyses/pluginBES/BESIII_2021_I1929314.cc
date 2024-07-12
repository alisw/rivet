// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> 4 charged particles (+pi0) cross sections
  class BESIII_2021_I1929314 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1929314);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<8;++ix) {
	book(_nCharged[ix], "TMP/nCharged_"+to_str(ix+1));
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
      if(ntotal==4) {
	if(nCount[211]==1 && nCount[-211]==1) {
	  if(nCount[321]==1 && nCount[-321]==1 )
	    _nCharged[0]->fill();
	  else if(nCount[2212]==1 && nCount[-2212]==1 )
	    _nCharged[3]->fill();
	}
	else if(nCount[321]==2 && nCount[-321]==2 )
	  _nCharged[1]->fill();
	else if (nCount[211]==2 && nCount[-211]==2 )
	  _nCharged[2]->fill();
      }
      else if(ntotal==5 && nCount[111]==1) {
	if(nCount[211]==1 && nCount[-211]==1) {
	  if(nCount[321]==1 && nCount[-321]==1 )
	    _nCharged[4]->fill();
	  else if(nCount[2212]==1 && nCount[-2212]==1 )
	    _nCharged[7]->fill();
	}
	else if(nCount[321]==2 && nCount[-321]==2 )
	  _nCharged[5]->fill();
	else if (nCount[211]==2 && nCount[-211]==2 )
	  _nCharged[6]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<8;++ix) {
	double sigma = _nCharged[ix]->val()*fact;
	double error = _nCharged[ix]->err()*fact;
	Scatter2D temphisto(refData(1, 1, ix+1));
	Scatter2DPtr  mult;
        book(mult, 1, 1, ix+1);
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
    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nCharged[8];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1929314);

}
