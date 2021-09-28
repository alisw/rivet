// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1980_I153511 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1980_I153511);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");

      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      
      // Thrust and sphericity
      declare(Sphericity(cfs), "Sphericity");
      // Book histograms
      unsigned int ihist=0;
      if      (fuzzyEquals(sqrtS()/GeV, 12., 1e-3)) {
	ihist=1;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 30., 1e-3)) {
	ihist=2;
      }
      else
	MSG_ERROR("Beam energy not supported!");

      book(_h_S ,   ihist, 1, 1);
      book(_h_A , 2+ihist, 1, 1);
      book(_h_x , 4+ihist, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h_S->fill(sphericity.sphericity());
      _h_A->fill(sphericity.aplanarity());
      for (const Particle& p : cfs.particles()) {
        const Vector3 mom3 = p.p3();
        const double mom = mom3.mod();
        const double xp = mom/meanBeamMom;
        _h_x->fill(xp);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_S);
      normalize(_h_A);
      normalize(_h_x);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_S,_h_A,_h_x;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1980_I153511);


}
