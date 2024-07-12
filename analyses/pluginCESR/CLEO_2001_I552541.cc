// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEO_2001_I552541 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2001_I552541);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_d_Dbar0[0], "/TMP/d_D0_low" );
      book(_d_Dbar0[1], "/TMP/d_D0_high");
      book(_d_Dm[0]   , "/TMP/d_Dm_low" );
      book(_d_Dm[1]   , "/TMP/d_Dm_high");
      book(_d_Lam[0]  , "/TMP/d_La_low" );
      book(_d_Lam[1]  , "/TMP/d_La_high");
      
      book(_n_Dbar0[0][0], "/TMP/d_D0_low_low"  );
      book(_n_Dbar0[0][1], "/TMP/d_D0_low_high" );
      book(_n_Dbar0[1][0], "/TMP/d_D0_high_low" );
      book(_n_Dbar0[1][1], "/TMP/d_D0_high_high");
      book(_n_Dm[0][0]   , "/TMP/d_Dm_low_low"  );
      book(_n_Dm[0][1]   , "/TMP/d_Dm_low_high" );
      book(_n_Dm[1][0]   , "/TMP/d_Dm_high_low" );
      book(_n_Dm[1][1]   , "/TMP/d_Dm_high_high");
      book(_n_Lam[0][0]  , "/TMP/d_La_low_low"  );
      book(_n_Lam[0][1]  , "/TMP/d_La_low_high" );
      book(_n_Lam[1][0]  , "/TMP/d_La_high_low" );
      book(_n_Lam[1][1]  , "/TMP/d_La_high_high");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles(Cuts::pid==-4122 or Cuts::pid==-411 or Cuts::pid==-421)) {
	long id1 = p.pid();
	double mom1 = p.p3().mod();
	if(mom1<2.3*GeV || mom1>5.*GeV) continue;
	bool high1 = mom1>3.3*GeV;
	if(id1==-4122) {
	  _d_Lam[high1]->fill();
	}
	else if(id1==-411) {
	  _d_Dm[high1]->fill();
	}
	else if(id1==-421) {
	  _d_Dbar0[high1]->fill();
	}
	for(const Particle & p2 : ufs.particles(Cuts::pid==4122)) {
	  if(p.p3().angle(p2.p3())<0.5*M_PI) continue;
	  double mom2 = p2.p3().mod();
	  if(mom2<2.3*GeV || mom2>5.*GeV) continue;
	  bool high2 = mom2>3.3*GeV;
	  if(id1==-4122) {
	    _n_Lam[high1][high2]->fill();
	  }
	  else if(id1==-411) {
	    _n_Dm[high1][high2]->fill();
	  }
	  else if(id1==-421) {
	    _n_Dbar0[high1][high2]->fill();
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter1D R_D0_low_low    = *_n_Dbar0[0][0]/ *_d_Dbar0[0];
      Scatter1D R_D0_low_high   = *_n_Dbar0[0][1]/ *_d_Dbar0[0];
      Scatter1D R_D0_high_low   = *_n_Dbar0[1][0]/ *_d_Dbar0[1];
      Scatter1D R_D0_high_high  = *_n_Dbar0[1][1]/ *_d_Dbar0[1];
      Scatter1D R_Dm_low_low    = *_n_Dm[0][0]   / *_d_Dm[0];
      Scatter1D R_Dm_low_high   = *_n_Dm[0][1]   / *_d_Dm[0];
      Scatter1D R_Dm_high_low   = *_n_Dm[1][0]   / *_d_Dm[1];
      Scatter1D R_Dm_high_high  = *_n_Dm[1][1]   / *_d_Dm[1];
      Scatter1D R_Lam_low_low   = *_n_Lam[0][0]  / *_d_Lam[0];
      Scatter1D R_Lam_low_high  = *_n_Lam[0][1]  / *_d_Lam[0];
      Scatter1D R_Lam_high_low  = *_n_Lam[1][0]  / *_d_Lam[1];
      Scatter1D R_Lam_high_high = *_n_Lam[1][1]  / *_d_Lam[1];
      for(unsigned int ix=3;ix<5;++ix) {
	for(unsigned int iy=1;iy<5;++iy) {
	  double num(0.),den(0.),num_err(0.),den_err(0.);
	  if(ix==3) {
	    if(iy==1) {
	      den     =  R_D0_low_low  .points()[0].x();
	      den_err =  R_D0_low_low  .points()[0].xErrAvg();
	    }
	    else if(iy==2) {
	      den     =  R_D0_high_low .points()[0].x();
	      den_err =  R_D0_high_low .points()[0].xErrAvg();
	    }
	    else if(iy==3) {
	      den     =  R_D0_low_low  .points()[0].x();
	      den_err =  R_D0_low_low  .points()[0].xErrAvg();
	    }
	    else if(iy==4) {
	      den     =  R_D0_high_high.points()[0].x();
	      den_err =  R_D0_high_high.points()[0].xErrAvg();
	    }
	  }
	  else if(ix==4) {
	    if(iy==1) {
	      den     =  R_Dm_low_low  .points()[0].x();
	      den_err =  R_Dm_low_low  .points()[0].xErrAvg();
	    }
	    else if(iy==2) {
	      den     =  R_Dm_high_low .points()[0].x();
	      den_err =  R_Dm_high_low .points()[0].xErrAvg();
	    }
	    else if(iy==3) {
	      den     =  R_Dm_low_low  .points()[0].x();
	      den_err =  R_Dm_low_low  .points()[0].xErrAvg();
	    }
	    else if(iy==4) {
	      den     =  R_Dm_high_high.points()[0].x();
	      den_err =  R_Dm_high_high.points()[0].xErrAvg();
	    }
	  }
	  if(iy==1) {
	    num     =  R_Lam_low_low  .points()[0].x();
	    num_err =  R_Lam_low_low  .points()[0].xErrAvg();
	  }
	  else if(iy==2) {
	    num     =  R_Lam_high_low .points()[0].x();
	    num_err =  R_Lam_high_low .points()[0].xErrAvg();
	  }
	  else if(iy==3) {
	    num     =  R_Lam_low_low  .points()[0].x();
	    num_err =  R_Lam_low_low  .points()[0].xErrAvg();
	  }
	  else if(iy==4) {
	    num     =  R_Lam_high_high.points()[0].x();
	    num_err =  R_Lam_high_high.points()[0].xErrAvg();
	  }
	  double val = num/den;
	  double err = val>=0. ? val*sqrt(sqr(num_err/num)+sqr(den_err/den)) : 0.;
	  Scatter2DPtr ratio;
	  book(ratio,ix, 1, iy);
	  Scatter2D temphisto(refData(ix, 1, iy));
	  const double x  = temphisto.point(0).x();
	  pair<double,double> ex = temphisto.point(0).xErrs();
	  ratio->addPoint(x, val, ex, make_pair(err,err));
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _d_Dbar0[2],_d_Dm[2],_d_Lam[2];
    CounterPtr _n_Dbar0[2][2],_n_Dm[2][2],_n_Lam[2][2];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_2001_I552541);


}
