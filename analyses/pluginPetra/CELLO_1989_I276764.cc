// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CELLO_1989_I276764 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CELLO_1989_I276764);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_gamma , 2, 1, 1);
      book(_h_pi0A  , 3, 1, 1);
      book(_h_pi0B  , 4, 1, 1);
      book(_h_eta   , 5, 1, 1);
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
      for (const Particle& p : apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==111 || Cuts::pid==221)) {
	double xE = p.E()/meanBeamMom;
	if(p.pid()==111) {
	  _h_pi0A->fill(xE);
	  _h_pi0B->fill(xE);
	}
	else
	  _h_eta->fill(xE);
      }
      for (const Particle& p : apply<FinalState>(event, "FS").particles(Cuts::pid==111 || Cuts::pid==22)) {
	double xE = p.E()/meanBeamMom;
	if(p.pid()==22)
	  _h_gamma->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact = sqr(sqrtS())/GeV2*crossSection()/microbarn/sumOfWeights();
      scale(_h_gamma, fact); 
      scale(_h_pi0A , fact); 
      scale(_h_pi0B , fact); 
      scale(_h_eta  , fact); 

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_gamma,_h_pi0A,_h_pi0B,_h_eta;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CELLO_1989_I276764);


}
