// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief average pT w.r.t. thrust and sphericity axes
  class PLUTO_1983_I191161 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1983_I191161);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      // Thrust
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      const Sphericity sphericity(fs);
      declare(sphericity, "Sphericity");
      _iBin=-1;
      if(fuzzyEquals(sqrtS(),7.7,1e-2))
	_iBin = 0;
      else if(fuzzyEquals(sqrtS(),9.4,1e-2))
	_iBin = 1;
      else if(fuzzyEquals(sqrtS(),12.,1e-2))
	_iBin = 2;
      else if(fuzzyEquals(sqrtS(),13.,1e-2))
	_iBin = 3;
      else if(fuzzyEquals(sqrtS(),17.,1e-2))
	_iBin = 4;
      else if(fuzzyEquals(sqrtS(),22.,1e-2))
	_iBin = 5;
      else if(fuzzyEquals(sqrtS(),27.6,1e-2))
	_iBin = 6;
      else if(fuzzyEquals(sqrtS(),30.8,1e-2))
	_iBin = 7;
      else
	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
      
      // Book histograms
      book(_p_thrust_pt     , 1, 1, 1);
      book(_p_thrust_pt2    , 1, 1, 2);
      book(_p_thrust_sum_pt , 1, 1, 3);
      book(_p_thrust_sum_pt2, 1, 1, 4);
      book(_p_sphere_pt     , 2, 1, 1);
      book(_p_sphere_pt2    , 2, 1, 2);
      book(_p_sphere_sum_pt , 2, 1, 3);
      book(_p_sphere_sum_pt2, 2, 1, 4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const FinalState & fs = apply<FinalState>(event, "FS");
      int nPart = fs.particles().size();
      double pT_T_sum(0.),pT2_T_sum(0.);
      double pT_S_sum(0.),pT2_S_sum(0.);
      int nCharged(0);
      for(const Particle & p : fs.particles()) {
        const Vector3 mom3 = p.p3();
        const double pTinT = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis()); 
        const double pT2_T = sqr(pTinT) + sqr(pToutT);
        const double pT2_S = sqr(pTinS) + sqr(pToutS);
	if(PID::isCharged(p.pid())) ++nCharged;
      	pT_T_sum  += sqrt(pT2_T);
      	pT2_T_sum +=      pT2_T ;
      	pT_S_sum  += sqrt(pT2_S);
      	pT2_S_sum +=      pT2_S ;
      }
      if(nCharged<4) vetoEvent;
      _p_thrust_pt      ->bins()[_iBin].fill(sqrtS(),pT_T_sum /nPart/MeV         );
      _p_thrust_pt2     ->bins()[_iBin].fill(sqrtS(),pT2_T_sum/nPart/1e3/sqr(MeV));
      _p_thrust_sum_pt  ->bins()[_iBin].fill(sqrtS(),pT_T_sum /GeV               );
      _p_thrust_sum_pt2 ->bins()[_iBin].fill(sqrtS(),pT2_T_sum/GeV               );
      _p_sphere_pt      ->bins()[_iBin].fill(sqrtS(),pT_S_sum /nPart/MeV         );
      _p_sphere_pt2     ->bins()[_iBin].fill(sqrtS(),pT2_S_sum/nPart/1e3/sqr(MeV));
      _p_sphere_sum_pt  ->bins()[_iBin].fill(sqrtS(),pT_S_sum /GeV               );
      _p_sphere_sum_pt2 ->bins()[_iBin].fill(sqrtS(),pT2_S_sum/GeV               );

    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    //@}


    /// @name Histograms
    //@{
    Profile1DPtr _p_thrust_pt, _p_thrust_pt2, _p_thrust_sum_pt, _p_thrust_sum_pt2;
    Profile1DPtr _p_sphere_pt, _p_sphere_pt2, _p_sphere_sum_pt, _p_sphere_sum_pt2;
    unsigned int _iBin;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1983_I191161);


}
