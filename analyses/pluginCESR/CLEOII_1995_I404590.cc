// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c+ spectrum
  class CLEOII_1995_I404590 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1995_I404590);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_x,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idXi = 4232;
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.468));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idXi)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x,1,false);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1995_I404590);

}
