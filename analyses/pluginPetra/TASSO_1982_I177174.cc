// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Charged particle spectra for a range of energies
  class TASSO_1982_I177174 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1982_I177174);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      if(isCompatibleWithSqrtS(12.)) {
	book(_h_x2,2,1,1);
	book(_h_x3,3,1,1);
      }
      else if (isCompatibleWithSqrtS(14.)) {
	book(_h_x1,1,1,1);
	book(_h_x2,2,1,2);
	book(_h_x3,3,1,2);
      }
      else if (isCompatibleWithSqrtS(22.)) {
	book(_h_x1,1,1,2);
	book(_h_x2,2,1,3);
	book(_h_x3,3,1,3);
      }
      else if (isCompatibleWithSqrtS(25.)) {
	book(_h_x2,2,1,4);
	book(_h_x3,3,1,4);
      }
      else if (isCompatibleWithSqrtS(30.)) {
	book(_h_x2,2,1,5);
	book(_h_x3,3,1,5);
      }
      else if (isCompatibleWithSqrtS(34.)) {
	book(_h_x2,2,1,6);
	book(_h_x3,3,1,6);
      }
      else if (isCompatibleWithSqrtS(35.)) {
        book(_h_x2, 2,1,7);
        book(_h_x3, 3,1,7);
      }

      if(inRange(sqrtS()/GeV, 29.9,36.7))
        book(_h_x1, 1,1,3);

      if(_h_x1==Histo1DPtr() && _h_x2==Histo1DPtr() && _h_x3==Histo1DPtr())
      	MSG_ERROR("Beam energy not supported!");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      for (const Particle& p : fs.particles()) {
	double xp = p.p3().mod()/meanBeamMom;
	if(_h_x1!=Histo1DPtr()) _h_x1->fill(xp);
	if(_h_x2!=Histo1DPtr()) _h_x2->fill(xp);
	if(_h_x3!=Histo1DPtr()) _h_x3->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_x1, crossSection()/microbarn/sumOfWeights()*sqr(sqrtS()));
      scale(_h_x2, crossSection()/microbarn/sumOfWeights()*sqr(sqrtS()));
      scale(_h_x3, 1./sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x1,_h_x2,_h_x3;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1982_I177174);


}
