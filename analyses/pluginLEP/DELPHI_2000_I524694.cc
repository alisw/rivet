// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DELPHI_2000_I524694 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2000_I524694);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      book(_histXpSigma,  1, 1, 1);
      book(_histXpLambda,  3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Get event weight for histo filling

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const int id = p.abspid();
        double xp = p.p3().mod()/meanBeamMom;
        switch (id) {
        case 3112:
          _histXpSigma->fill(xp);
	  break;
        case 3124:
          _histXpLambda->fill(xp);
	  break;
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = 1./sumOfWeights();
      scale(_histXpSigma , fact); 
      scale(_histXpLambda, fact); 
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histXpSigma;
    Histo1DPtr _histXpLambda;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_2000_I524694);


}
