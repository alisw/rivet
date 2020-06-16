// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* production
  class OPAL_1995_I382219 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_1995_I382219);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      book(_h_Xe_Ds  , 3, 1, 1);
      book(_h_Xe_Ds_b, 4, 1, 1);
      book(_h_Xe_Ds_c, 5, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0/GeV;
      // check if b hadrons or not
      unsigned int nB= (filter_select(ufs.particles(), isBottomHadron)).size();
      // Accept all D*+- decays. 
      for (const Particle& p : filter_select(ufs.particles(), Cuts::abspid==PID::DSTARPLUS)) {
	// Scaled energy.
	const double energy = p.E()/GeV;
	const double scaledEnergy = energy/meanBeamMom;
	_h_Xe_Ds->fill(scaledEnergy);
	if(nB==0)
	  _h_Xe_Ds_c->fill(scaledEnergy);
	else
	  _h_Xe_Ds_b->fill(scaledEnergy);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // brs for D*+/- -> D0 pi+/- and D0->K+pi-
      double br = 0.677*0.03950;
      scale(_h_Xe_Ds  , 1./sumOfWeights()*br);
      scale(_h_Xe_Ds_b, 1./sumOfWeights()*br);
      scale(_h_Xe_Ds_c, 1./sumOfWeights()*br);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Xe_Ds,_h_Xe_Ds_b,_h_Xe_Ds_c;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1995_I382219);


}
