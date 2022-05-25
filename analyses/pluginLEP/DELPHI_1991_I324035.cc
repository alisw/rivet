// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged particle multiplicities in different regions
  class DELPHI_1991_I324035 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1991_I324035);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");

      // Book histograms
      book(_h_all_05  ,  1, 1, 1);
      book(_h_all_10  ,  2, 1, 1);
      book(_h_all_15  ,  3, 1, 1);
      book(_h_all_20  ,  4, 1, 1);
      book(_h_all_all ,  5, 1, 1);
      book(_h_hemi_05 ,  6, 1, 1);
      book(_h_hemi_10 ,  7, 1, 1);
      book(_h_hemi_15 ,  8, 1, 1);
      book(_h_hemi_20 ,  9, 1, 1);
      book(_h_hemi_30 , 10, 1, 1);
      book(_h_hemi_40 , 11, 1, 1);
      book(_h_hemi_50 , 12, 1, 1);
      book(_h_hemi_all, 13, 1, 1);
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

      // Thrusts
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis=thrust.thrustAxis();

      unsigned int n_all_05(0),n_all_10(0),n_all_15(0),n_all_20(0),n_all_all(0);
      unsigned int n_pos_05(0),n_pos_10(0),n_pos_15(0),n_pos_20(0),n_pos_30(0),n_pos_40(0),n_pos_50(0),n_pos_all(0);
      unsigned int n_neg_05(0),n_neg_10(0),n_neg_15(0),n_neg_20(0),n_neg_30(0),n_neg_40(0),n_neg_50(0),n_neg_all(0);
      for(const Particle &p : fs.particles()) {
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
        const double momT = dot(axis, mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
	++n_all_all;
	if(abs(rapidityT)<0.5) {
	  ++n_all_05;
	  if(rapidityT>0)
	    ++n_pos_05;
	  else
	    ++n_neg_05;
	}
	if(abs(rapidityT)<1.0) {
	  ++n_all_10;
	  if(rapidityT>0)
	    ++n_pos_10;
	  else
	    ++n_neg_10;
	}
	if(abs(rapidityT)<1.5) {
	  ++n_all_15;
	  if(rapidityT>0)
	    ++n_pos_15;
	  else
	    ++n_neg_15;
	}
	if(abs(rapidityT)<2.0) {
	  ++n_all_20;
	  if(rapidityT>0)
	    ++n_pos_20;
	  else
	    ++n_neg_20;
	}
	if(abs(rapidityT)<3.0) {
	  if(rapidityT>0)
	    ++n_pos_30;
	  else
	    ++n_neg_30;
	}
	if(abs(rapidityT)<4.0) {
	  if(rapidityT>0)
	    ++n_pos_40;
	  else
	    ++n_neg_40;
	}
	if(abs(rapidityT)<5.0) {
	  if(rapidityT>0)
	    ++n_pos_50;
	  else
	    ++n_neg_50;
	}
	if(rapidityT>0)
	  ++n_pos_all;
	else
	  ++n_neg_all;
      }
      _h_all_05 ->fill(n_all_05 );
      _h_all_10 ->fill(n_all_10 );
      _h_all_15 ->fill(n_all_15 );
      _h_all_20 ->fill(n_all_20 );
      _h_all_all->fill(n_all_all);
      _h_hemi_05 ->fill(n_pos_05 );
      _h_hemi_10 ->fill(n_pos_10 );
      _h_hemi_15 ->fill(n_pos_15 );
      _h_hemi_20 ->fill(n_pos_20 );
      _h_hemi_30 ->fill(n_pos_30 );
      _h_hemi_40 ->fill(n_pos_40 );
      _h_hemi_50 ->fill(n_pos_50 );
      _h_hemi_all->fill(n_pos_all);
      _h_hemi_05 ->fill(n_neg_05 );
      _h_hemi_10 ->fill(n_neg_10 );
      _h_hemi_15 ->fill(n_neg_15 );
      _h_hemi_20 ->fill(n_neg_20 );
      _h_hemi_30 ->fill(n_neg_30 );
      _h_hemi_40 ->fill(n_neg_40 );
      _h_hemi_50 ->fill(n_neg_50 );
      _h_hemi_all->fill(n_neg_all);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize( _h_all_05  , 1000.);
      normalize( _h_all_10  , 1000.);
      normalize( _h_all_15  , 1000.);
      normalize( _h_all_20  , 1000.);
      normalize( _h_all_all , 2000.);
      normalize( _h_hemi_05 , 1000.);
      normalize( _h_hemi_10 , 1000.);
      normalize( _h_hemi_15 , 1000.);
      normalize( _h_hemi_20 , 1000.);
      normalize( _h_hemi_30 , 1000.);
      normalize( _h_hemi_40 , 1000.);
      normalize( _h_hemi_50 , 1000.);
      normalize( _h_hemi_all, 1000.);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_all_05, _h_all_10, _h_all_15, _h_all_20, _h_all_all;
    Histo1DPtr _h_hemi_05, _h_hemi_10, _h_hemi_15, _h_hemi_20,
      _h_hemi_30, _h_hemi_40, _h_hemi_50, _h_hemi_all;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_1991_I324035);


}
