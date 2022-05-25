// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Omega_c*0 production
  class BABAR_2006_I724557 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I724557);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histos
      book(_c_Omega     ,"TMP/Omega"    ,1,0.,1.);
      book(_c_Omega_Star,"TMP/OmegaStar",1,0.,1.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & omega : apply<UnstableParticles>(event,"UFS").particles(Cuts::abspid==4332 or
										   Cuts::abspid==4334)) {
	double xp = omega.momentum().p3().mod()/sqrt(0.25*sqr(sqrtS())-sqr(omega.mass()));
	if(xp<0.5) continue;
	if(omega.abspid()==4332)
	  _c_Omega->fill(0.5);
	else
	  _c_Omega_Star->fill(0.5);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter2DPtr ratio;
      book(ratio,1,1,1);
      divide(_c_Omega_Star,_c_Omega,ratio);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _c_Omega,_c_Omega_Star;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I724557);

}
