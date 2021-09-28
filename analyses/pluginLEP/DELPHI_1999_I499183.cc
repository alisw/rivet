// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief event shapes at 133, 161 172 and 183 GeV
  class DELPHI_1999_I499183 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1999_I499183);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      unsigned int offset = 0;
      int offset2 = 0;

      if (fuzzyEquals(sqrtS()/GeV, 133  , 1E-3)) {
	offset  = 0;			   
	offset2 = 1;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 161  , 1E-3)) {
	offset  = 0;			   
	offset2 = 2;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 172  , 1E-3)) {
	offset  = 0;			   
	offset2 = 3;			   
      }					   
      else if (fuzzyEquals(sqrtS()/GeV, 183  , 1E-3)) {
	offset  = 1;			   
	offset2 = 1;			   
      }
      
      book(_h_thrust          , 13+offset, 1, offset2);
      book(_h_major           , 15+offset, 1, offset2);
      book(_h_minor           , 17+offset, 1, offset2);
      book(_h_oblateness      , 19+offset, 1, offset2);
      book(_h_sphericity      , 21+offset, 1, offset2);
      book(_h_planarity       , 23+offset, 1, offset2);
      book(_h_aplanarity      , 25+offset, 1, offset2);
      book(_h_heavy_jet_mass  , 27+offset, 1, offset2);
      book(_h_light_jet_mass  , 29+offset, 1, offset2);
      book(_h_diff_jet_mass   , 31+offset, 1, offset2);
      book(_h_wide_broading   , 33+offset, 1, offset2);
      book(_h_narrow_broading , 35+offset, 1, offset2);
      book(_h_total_broading  , 37+offset, 1, offset2);
      book(_h_diff_broading   , 39+offset, 1, offset2);
      book(_h_CParam          , 41+offset, 1, offset2);
      book(_h_DParam          , 43+offset, 1, offset2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // thrust related observables
      _h_thrust    ->fill(1.-thrust.thrust()  );
      _h_major     ->fill(thrust.thrustMajor());
      _h_minor     ->fill(thrust.thrustMinor());
      _h_oblateness->fill(thrust.oblateness() );

      // sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h_sphericity->fill(sphericity.sphericity());
      _h_planarity ->fill(sphericity.planarity() );
      _h_aplanarity->fill(sphericity.aplanarity());
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      // standard jet masses
      _h_heavy_jet_mass->fill(hemi.scaledM2high());
      _h_light_jet_mass->fill(hemi.scaledM2low() );
      _h_diff_jet_mass ->fill(hemi.scaledM2diff());
      // jet broadening
      _h_wide_broading  ->fill(hemi.Bmax() );
      _h_narrow_broading->fill(hemi.Bmin() );
      _h_total_broading ->fill(hemi.Bsum() );
      _h_diff_broading  ->fill(hemi.Bdiff());
      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h_CParam->fill(parisi.C());
      _h_DParam->fill(parisi.D());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_thrust          );
      normalize(_h_major           );
      normalize(_h_minor           );
      normalize(_h_sphericity      );
      normalize(_h_planarity       );
      normalize(_h_aplanarity       );
      normalize(_h_oblateness      );
      normalize(_h_heavy_jet_mass  );
      normalize(_h_light_jet_mass  );
      normalize(_h_diff_jet_mass   );
      normalize(_h_wide_broading   );
      normalize(_h_narrow_broading );
      normalize(_h_total_broading  );
      normalize(_h_diff_broading   );
      normalize(_h_CParam   );
      normalize(_h_DParam   );

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust,_h_major,_h_minor;
    Histo1DPtr _h_sphericity,_h_planarity,_h_aplanarity,_h_oblateness;
    Histo1DPtr _h_heavy_jet_mass,_h_light_jet_mass,_h_diff_jet_mass;
    Histo1DPtr _h_wide_broading,_h_narrow_broading,_h_total_broading,_h_diff_broading;
    Histo1DPtr _h_CParam,_h_DParam;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1999_I499183);


}
