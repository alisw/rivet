// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c' spectrum
  class CLEOII_1999_I478217 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1999_I478217);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_Xi_c,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.578));
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==4312 or Cuts::abspid==4322)) {
	double xp = p.momentum().p3().mod()/Pmax;
        _h_Xi_c->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Xi_c);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Xi_c;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1999_I478217);

}
