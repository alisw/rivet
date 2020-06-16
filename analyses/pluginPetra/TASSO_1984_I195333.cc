// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1984_I195333 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1984_I195333);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      declare(Sphericity(fs), "Sphericity");
      declare(Thrust(fs), "Thrust");
      
      // counters for R
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_h_weight, "/TMP/HWeight");
      unsigned int iloc(0);
      if(fuzzyEquals(sqrtS()/GeV, 14 , 1E-3))
	iloc = 1;
      else if(fuzzyEquals(sqrtS()/GeV, 22 , 1E-3))
	iloc = 2;
      else if(fuzzyEquals(sqrtS()/GeV, 34 , 1E-3))
	iloc = 3;
      if(iloc!=0) {
	book(_h_mult,  3,1,iloc);
	book(_h_p,  5,1,iloc);
	book(_h_xp,  6,1,iloc);
	book(_h_pl,  7,1,iloc);
	book(_h_pt,  8,1,iloc);
	book(_h_pt2,  9,1,iloc);
	book(_h_xl, 10,1,iloc);
	book(_h_xT, 11,1,iloc);
	book(_h_S, 12,1,iloc);
	book(_h_T, 13,1,iloc);
	book(_h_y, 14,1,iloc);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      unsigned int nCharged(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(PID::isCharged(p.pid())) ++nCharged;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22]) {	
	_c_muons->fill();
	return;
      }
      // everything else
      _c_hadrons->fill();
      _h_weight->fill();
      _n_charged.fill(nCharged);
      _n_total.fill(ntotal);
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _thrust.fill(thrust.thrust());
      // sphericity
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _sphericity.fill(sphericity.sphericity());
      // global distributions
      if(_h_mult) _h_mult->fill(nCharged);
      if(_h_S)    _h_S   ->fill(sphericity.sphericity());
      if(_h_T)    _h_T   ->fill(thrust.thrust());
      // single particle distributions
      for (const Particle& p : fs.particles()) {
	if(!PID::isCharged(p.pid())) continue;
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	_p_total.fill(pp);
	if(_h_p)  _h_p ->fill(pp);
	if(_h_xp) _h_xp->fill(2.*pp/sqrtS());
        const double mom = dot(sphericity.sphericityAxis(), mom3);
	_p_l.fill(fabs(mom));
	if(_h_pl) _h_pl->fill(fabs(mom));
	if(_h_xl) _h_xl->fill(2.*fabs(mom)/sqrtS());
        const double pTin = dot(mom3, sphericity.sphericityMajorAxis());
	_pt2_in.fill(sqr(pTin));
        const double pTout = dot(mom3, sphericity.sphericityMinorAxis());
	_pt2_out.fill(sqr(pTout));
        double pT = sqr(pTin) + sqr(pTout);
	_pt2.fill(pT);
	if(_h_pt2) _h_pt2->fill(pT);
	pT=sqrt(pT);
	_pt.fill(pT);
	if(_h_pt) _h_pt->fill(pT);
	if(_h_xT) _h_xT->fill(2.*pT/sqrtS());
	if(_h_y) {
	  const double rap = 0.5 * log((p.E() + mom) /
				       (p.E() - mom));
	  _h_y->fill(fabs(rap));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter1D R = *_c_hadrons/ *_c_muons;
      double              rval = R.point(0).x();
      pair<double,double> rerr = R.point(0).xErrs();
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr hadrons;
      book(hadrons, "sigma_hadrons");
      Scatter2DPtr muons;
      book(muons, "sigma_muons"  );
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, rval, ex, rerr);
	  hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
	  muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	  hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	  muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
      // charged particle multiplicity distribution
      if(_h_mult) normalize(_h_mult,2.);
      for(unsigned int iy=1;iy<12;++iy) {
	double value = 0.0, error = 0.0;
	if(iy==1) {
	  value = _n_charged.xMean();
	  error = _n_charged.xStdErr();
	}
	else if(iy==2) {
	  double num = _n_charged.xMean();
	  double den =   _n_total.xMean();
	  value = num/den;
	  error = value*sqrt(sqr(_n_charged.xStdErr()/num)+sqr(_n_total.xStdErr()/den));
	}
	else if(iy==3) {
	  value = _n_charged.xStdDev();
	  error = _n_charged.xStdErr();
	}
	else if(iy==4) {
	  value = _sphericity.xMean();
	  error = _sphericity.xStdErr();
	}
	else if(iy==5) {
	  value = _thrust.xMean();
	  error = _thrust.xStdErr();
	}
	else if(iy==6) {
	  value = _p_total.xMean();
	  error = _p_total.xStdErr();
	}
	else if(iy==7) {
	  value = _p_l.xMean();
	  error = _p_l.xStdErr();
	}
	else if(iy==8) {
	  value = _pt.xMean();
	  error = _pt.xStdErr();
	}
	else if(iy==9) {
	  value = _pt2.xMean();
	  error = _pt2.xStdErr();
	}
	else if(iy==10) {
	  value = _pt2_in.xMean();
	  error = _pt2_in.xStdErr();
	}
	else if(iy==11) {
	  value = _pt2_out.xMean();
	  error = _pt2_out.xStdErr();
	}
	Scatter2D temphisto(refData(4, 1, iy));
	Scatter2DPtr mult;
	book(mult, 4, 1, iy);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult   ->addPoint(x, value, ex, make_pair(error,error));
	  }
	  else {
	    mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
      // scale the distributions
      scale(_h_p  ,1./_h_weight->sumW());
      scale(_h_xp ,1./_h_weight->sumW());
      scale(_h_pl ,1./_h_weight->sumW());
      scale(_h_pt ,1./_h_weight->sumW());
      scale(_h_pt2,1./_h_weight->sumW());
      scale(_h_xl ,1./_h_weight->sumW());
      scale(_h_xT ,1./_h_weight->sumW());
      scale(_h_S  ,1./_h_weight->sumW());
      scale(_h_T  ,1./_h_weight->sumW());
      scale(_h_y  ,1./_h_weight->sumW());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_mult,_h_p,_h_xp,_h_pl,_h_pt,_h_pt2,_h_xl,_h_xT,_h_S,_h_T,_h_y;
    CounterPtr _c_hadrons, _c_muons;
    YODA::Dbn1D _n_charged,_n_total,_sphericity,_thrust,_p_total,
      _p_l,_pt,_pt2,_pt2_in,_pt2_out;
    CounterPtr  _h_weight;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1984_I195333);

}
