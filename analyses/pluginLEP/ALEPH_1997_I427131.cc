// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi0 spectrum
  class ALEPH_1997_I427131 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1997_I427131);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // // Initialise and register projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      declare(Sphericity(fs), "Sphericity");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_x  , 2, 1, 2);
      book(_h_in , 4, 1, 1);
      book(_h_out, 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.particles().size() < 2) {
      	MSG_DEBUG("Failed ncharged cut");
      	vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==111) ) {
	const Vector3 mom3 = p.p3();
	const double pTinS  = abs(dot(mom3, sphericity.sphericityMajorAxis()));
	const double pToutS = abs(dot(mom3, sphericity.sphericityMinorAxis()));
       	_h_x  ->fill(mom3.mod()/meanBeamMom);
       	_h_in ->fill(pTinS                 );
       	_h_out->fill(pToutS                );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_x  , 1./sumOfWeights());
      scale(_h_in , 1./sumOfWeights());
      scale(_h_out, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x, _h_in, _h_out;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALEPH_1997_I427131);


}
