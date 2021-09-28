// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TPC_1984_I200105 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TPC_1984_I200105);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      const FinalState cfs(Cuts::charge != 0);
      declare(cfs, "FS");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_x , 1, 1, 1);
      book(_h_pt, 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // require 5 charged particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.size()<3) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // calc thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis = thrust.thrustAxis();
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==333)) {
        double xE = p.E()/meanBeamMom;
	double pT2 = p.p3().mod2()-sqr(axis.dot(p.p3()));
	_h_x ->fill(xE );
	_h_pt->fill(pT2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_x, 1./sumOfWeights());
      scale(_h_pt, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x,_h_pt;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TPC_1984_I200105);


}
