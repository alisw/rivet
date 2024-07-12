// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charm hadrons at 91 GeV
  class DELPHI_1993_I356732 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1993_I356732);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");

      book(_h_Xe_Ds  , 1, 1, 1);
      book(_h_Xe_D0  , 3, 1, 1);
      book(_h_Xe_Dp  , 5, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      
      // Accept all D*+- decays. 
      for (const Particle& p : ufs.particles(Cuts::abspid==PID::DSTARPLUS ||
					     Cuts::abspid==PID::DPLUS ||
					     Cuts::abspid==PID::D0)) {
	// Scaled energy.
	const double energy = p.E()/GeV;
	const double scaledEnergy = 2.*energy/sqrtS();
	if(p.abspid()==PID::DSTARPLUS)
	   _h_Xe_Ds->fill(scaledEnergy);
	else if(p.abspid()==PID::DPLUS)
	   _h_Xe_Dp->fill(scaledEnergy);
	else
	   _h_Xe_D0->fill(scaledEnergy);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Xe_Ds  , 1./sumOfWeights());
      scale(_h_Xe_Dp  , 1./sumOfWeights());
      scale(_h_Xe_D0  , 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Xe_Ds,_h_Xe_D0,_h_Xe_Dp;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_1993_I356732);

}
