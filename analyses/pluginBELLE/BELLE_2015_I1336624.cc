// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2015_I1336624 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2015_I1336624);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_1S     , "/TMP/1S");
      book(_c_2S     , "/TMP/2S");
      book(_c_3S     , "/TMP/3S");
      book(_c_muons, "/TMP/sigma_muons");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // analyse the final state
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // intermediates
      bool isBottom(false);
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	// check for bottom hadrons
        if (PID::isBottomHadron(p.pid())) {
	  isBottom = true;
	  break;
	}
	// upsilon + pi+pi-
	if(p.children().empty()) continue;
	if(p.pid() !=   553  &&
	   p.pid() != 100553 &&
	   p.pid() != 200553 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==211) {
	    continue;
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched) {
	  if(nRes[211]==1 && nRes[-211]==1 ) {
	    if(p.pid()==553)
	      _c_1S->fill();
	    if(p.pid()==100553)
	      _c_2S->fill();
	    if(p.pid()==200553)
	      _c_3S->fill();
	  }
	}
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // open bottom
      else if(isBottom) {
	_c_hadrons->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // loop over histos to be filled
      for(unsigned int ix=1;ix<3;++ix) {
	Scatter1D R;
	for(unsigned int iy=1;iy<4;++iy) {
	  if(ix==2 && iy!=1) continue;
	  if(ix==1) {
	    if(iy==1) {
	      R = *_c_1S/ *_c_muons;
	    }
	    else if(iy==2) {
	      R = *_c_2S/ *_c_muons;
	    }
	    else {
	      R = *_c_3S/ *_c_muons;
	    }
	  }
	  else if(ix==2) {
	    R = *_c_hadrons/ *_c_muons;
	  }
	  double              rval = R.point(0).x();
	  pair<double,double> rerr = R.point(0).xErrs();
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr     mult;
	  book(mult, ix, 1, iy);
	  for (size_t b = 0; b < temphisto.numPoints(); b++) {
	    const double x  = temphisto.point(b).x();
	    pair<double,double> ex = temphisto.point(b).xErrs();
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	      mult   ->addPoint(x, rval, ex, rerr);
	    }
	    else {
	      mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_hadrons, _c_muons, _c_1S, _c_2S, _c_3S;
    //@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2015_I1336624);

}
