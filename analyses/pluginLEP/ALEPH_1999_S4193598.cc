// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  class ALEPH_1999_S4193598 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ALEPH_1999_S4193598()
      : Analysis("ALEPH_1999_S4193598")
    { }

    //@}


  public:

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(ChargedFinalState(), "CFS");

      book(_h_Xe_Ds ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
 
      // Trigger condition
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() < 5) vetoEvent;

      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0/GeV;

      // Accept all D*+- decays.
      for (const Particle& p : filter_select(ufs.particles(), Cuts::abspid==PID::DSTARPLUS)) {
          // Scaled energy.
          const double energy = p.E()/GeV;
          const double scaledEnergy = energy/meanBeamMom;
          _h_Xe_Ds->fill(scaledEnergy);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // brs for D*+/- -> D0 pi+/- and D0->K+pi-
      double br = 0.677*0.03950;
      scale(_h_Xe_Ds, 1./sumOfWeights()*br*1000.);
    }

  private:

    Histo1DPtr _h_Xe_Ds;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_1999_S4193598);

}
