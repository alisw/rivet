// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Lambda production at 29 GeV
  class HRS_1986_I17781 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1986_I17781);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");
      book(_h_spect,1,1,1);
      book(_h_rap  ,2,1,1);
      book(_h_mult ,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      // Thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==3122)) {
      	double xE = 2.*p.E()/sqrtS();
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
	_h_spect->fill(xE,1./beta);
	_h_rap->fill(abs(rapidityT));
	_h_mult->fill(sqrtS());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _h_spect, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale( _h_rap  , 1./sumOfWeights());
      scale( _h_mult , 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_spect,_h_rap,_h_mult;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1986_I17781);

}
