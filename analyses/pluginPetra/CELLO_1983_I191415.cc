// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi0 and gamma spectra at 14, 22 and 34
  class CELLO_1983_I191415 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1983_I191415);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      unsigned int iloc(0);
      if(isCompatibleWithSqrtS(14.)) {
	iloc=1;
      }
      else if (isCompatibleWithSqrtS(22.)) {
	iloc=2;
      }
      else if (isCompatibleWithSqrtS(34.)) {
	iloc=3;
      }
      else
	MSG_ERROR("Beam energy not supported!");
      // Book histograms
      book(_h_gamma, iloc  , 1, 1);
      book(_h_pi0  , iloc+3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // at least 5 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      if (numParticles < 5) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state to get particle spectra
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==111)) {
	double xE = p.E()/meanBeamMom;
	_h_pi0->fill(xE);
      }
      for (const Particle& p : apply<FinalState>(event, "FS").particles(Cuts::pid==22)) {
	double xE = p.E()/meanBeamMom;
	_h_gamma->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact = sqr(sqrtS())/GeV2*crossSection()/microbarn/sumOfWeights();
      scale(_h_gamma, fact); 
      scale(_h_pi0  , fact); 

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_gamma,_h_pi0;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CELLO_1983_I191415);


}
