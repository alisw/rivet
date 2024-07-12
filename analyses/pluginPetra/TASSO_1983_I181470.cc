// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, K and proton spectra at 14,22 and 34 GeV
  class TASSO_1983_I181470 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1983_I181470);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");

      vector<int> hist1,hist2;
      sqs = 1.;
      if(isCompatibleWithSqrtS(14.)) {
	hist1 = {19,21,23};
	hist2 = {20,22,24};
	sqs = 14.;
      }
      else if (isCompatibleWithSqrtS(22.)) {
	hist1 = {25,27,11};
	hist2 = {26,10,12};
	sqs = 22.;
      }
      else if (isCompatibleWithSqrtS(34.)) {
	hist1 = {13,15,17};
	hist2 = {14,16,18};
	sqs = 34.;
      }
      else
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
          << " doesn't match any available analysis energy .");
      
      book(_h_p_pi , hist1[0],1,1);
      book(_h_p_K  , hist1[1],1,1);
      book(_h_p_p  , hist1[2],1,1);
      book(_h_x_pi , hist2[0],1,1);
      book(_h_x_K  , hist2[1],1,1);
      book(_h_x_p  , hist2[2],1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
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
	double xE = p.E()/meanBeamMom;
	if(abs(p.pid())==211) {
	  _h_p_pi->fill(p.p3().mod());
	  _h_x_pi->fill(xE          );
	}
	else if(abs(p.pid())==321) {
	  _h_p_K->fill(p.p3().mod());
	  _h_x_K->fill(xE          );
	}
	else if(abs(p.pid())==2212) {
	  _h_p_p->fill(p.p3().mod());
	  _h_x_p->fill(xE          );
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact1 = crossSection()/nanobarn/sumOfWeights();
      double fact2 = sqr(sqs)/GeV2*crossSection()/microbarn/sumOfWeights();
      
      scale(_h_p_pi, fact1); 
      scale(_h_p_K , fact1); 
      scale(_h_p_p , fact1);
      
      scale(_h_x_pi, fact2); 
      scale(_h_x_K , fact2); 
      scale(_h_x_p , fact2); 
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr   _h_p_pi,_h_p_K,_h_p_p;
    Histo1DPtr   _h_x_pi,_h_x_K,_h_x_p;
    double sqs;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1983_I181470);

}
