// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charm meson spectra at 29 GeV
  class HRS_1988_I23360 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1988_I23360);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      // Book histograms
      book(_h_Dstar1,1,1,1);
      book(_h_Dstar2,1,1,2);
      book(_h_D0    ,2,1,1);
      book(_h_Dp    ,2,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      int nch = cfs.particles().size();
      if(nch<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==413 || Cuts::abspid==421 || Cuts::abspid==411 )) {
         double xE = p.E()/meanBeamMom;
	 int id = abs(p.pid());
	 if(id==413) {
	   _h_Dstar1->fill(xE);
	   _h_Dstar2->fill(xE);
	 }
	 else if(id==421) {
	   _h_D0->fill(xE);
	 }
	 else if(id==411) {
	   _h_Dp->fill(xE);
	 }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = sqr(sqrtS())*crossSection()/sumOfWeights()/microbarn;
      scale(_h_Dstar1, fact);
      normalize(_h_Dstar2);
      scale(_h_D0,fact);
      scale(_h_Dp,fact);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Dstar1,_h_Dstar2,_h_D0,_h_Dp;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(HRS_1988_I23360);


}
