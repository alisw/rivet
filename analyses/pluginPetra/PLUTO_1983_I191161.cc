// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief average pT w.r.t. thrust and sphericity axes
  class PLUTO_1983_I191161 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1983_I191161);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      declare(ChargedFinalState(), "CFS");
      // Thrust
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      const Sphericity sphericity(fs);
      declare(sphericity, "Sphericity");
      sqs = 1.0;
      if(isCompatibleWithSqrtS(7.7))
	sqs = 7.7;
      else if(isCompatibleWithSqrtS(9.4))
	sqs = 9.4;
      else if(isCompatibleWithSqrtS(12.))
	sqs = 12.;
      else if(isCompatibleWithSqrtS(13.))
	sqs = 13.;
      else if(isCompatibleWithSqrtS(17.))
	sqs = 17.;
      else if(isCompatibleWithSqrtS(22.))
	sqs = 22.;
      else if(isCompatibleWithSqrtS(27.6))
	sqs = 27.6;
      else if(isCompatibleWithSqrtS(30.8))
	sqs = 30.8;
      else
	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // at least 4 charged particles
      if(apply<ChargedFinalState>(event, "CFS").particles().size()<4) vetoEvent;
      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const FinalState & fs = apply<FinalState>(event, "FS");
      // remove pi0->gammagamma decay products and replace with the pi0s
      // needed to get the average pTs right
      Particles fsParticles;
      set<ConstGenParticlePtr> pi0;
      for(const Particle & p : fs.particles()) {
	if((p.pid()!=PID::PHOTON && p.abspid()!=PID::ELECTRON)|| p.parents().empty() || p.parents()[0].pid()!=PID::PI0)
	  fsParticles.push_back(p);
	else {
	  if(pi0.find(p.parents()[0].genParticle())==pi0.end())
	    fsParticles.push_back(p.parents()[0]);
	}
      }
      double pT_T_sum(0.),pT2_T_sum(0.);
      double pT_S_sum(0.),pT2_S_sum(0.);
      for(const Particle & p : fsParticles) {
        const Vector3 mom3 = p.p3();
        const double pTinT = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis()); 
        const double pT2_T = sqr(pTinT) + sqr(pToutT);
        const double pT2_S = sqr(pTinS) + sqr(pToutS);
        const double pT_T  = sqrt(pT2_T);
        const double pT_S  = sqrt(pT2_S);
      	pT_T_sum  += sqrt(pT2_T);
      	pT2_T_sum +=      pT2_T ;
      	pT_S_sum  += sqrt(pT2_S);
      	pT2_S_sum +=      pT2_S ;
	_p_thrust_pt      .fill(pT_T /MeV         );
	_p_thrust_pt2     .fill(pT2_T/1e3/sqr(MeV));
	_p_sphere_pt      .fill(pT_S /MeV         );
	_p_sphere_pt2     .fill(pT2_S/1e3/sqr(MeV));
      }
      _p_thrust_sum_pt  .fill(pT_T_sum /GeV               );
      _p_thrust_sum_pt2 .fill(pT2_T_sum/GeV               );
      _p_sphere_sum_pt  .fill(pT_S_sum /GeV               );
      _p_sphere_sum_pt2 .fill(pT2_S_sum/GeV               );
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
	for(unsigned int iy=1;iy<5;++iy) {
	  double value = 0.0, error = 0.0;
	  if(ix==1) {
	    if(iy==1) {
	      value = _p_thrust_pt.xMean();
	      error = _p_thrust_pt.xStdErr();
	    }
	    else if(iy==2) {
	      value = _p_thrust_pt2.xMean();
	      error = _p_thrust_pt2.xStdErr();
	    }
	    else if(iy==3) {
	      value = _p_thrust_sum_pt.xMean();
	      error = _p_thrust_sum_pt.xStdErr();
	    } 
	    else if(iy==4) {
	      value = _p_thrust_sum_pt2.xMean();
	      error = _p_thrust_sum_pt2.xStdErr();
	    } 
	  }
	  else {
	    if(iy==1) {
	      value = _p_sphere_pt.xMean();
	      error = _p_sphere_pt.xStdErr();
	    }
	    else if(iy==2) {
	      value = _p_sphere_pt2.xMean();
	      error = _p_sphere_pt2.xStdErr();
	    }
	    else if(iy==3) {
	      value = _p_sphere_sum_pt.xMean();
	      error = _p_sphere_sum_pt.xStdErr();
	    } 
	    else if(iy==4) {
	      value = _p_sphere_sum_pt2.xMean();
	      error = _p_sphere_sum_pt2.xStdErr();
	    } 
	  }
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr mult;
	  book(mult, ix, 1, iy);
	  for (size_t b = 0; b < temphisto.numPoints(); b++) {
	    const double x  = temphisto.point(b).x();
	    pair<double,double> ex = temphisto.point(b).xErrs();
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    if (inRange(sqs, x-ex2.first, x+ex2.second)) {
	      mult   ->addPoint(x, value, ex, make_pair(error,error));
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
    YODA::Dbn1D _p_thrust_pt, _p_thrust_pt2, _p_thrust_sum_pt, _p_thrust_sum_pt2;
    YODA::Dbn1D _p_sphere_pt, _p_sphere_pt2, _p_sphere_sum_pt, _p_sphere_sum_pt2;
    double sqs;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(PLUTO_1983_I191161);


}
