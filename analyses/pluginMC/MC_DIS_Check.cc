// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief A simple analysis to illustrate how to use the
  /// DISKinematics projection together with different options.
  class MC_DIS_Check : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_DIS_Check);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections. Note that the definition
      // of the scattered lepton can be influenced by sepcifying
      // options as declared in the .info file.
      DISLepton lepton(options());
      declare(lepton, "Lepton");
      declare(DISKinematics(lepton), "Kinematics");

      // Book histograms
	
      _hist_Q2  = bookHisto1D("Q2",logspace(100,0.1, 1000.0));
      _hist_y   = bookHisto1D("y",100,0.,1.);
      _hist_x   = bookHisto1D("xBj",logspace(100,0.00001, 1.0));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get the DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) return;
      double x  = dk.x();
      double y  = dk.y();
      double Q2 = dk.Q2();
      // Weight of the event
      const double weight = event.weight();
      _hist_Q2->fill(Q2,weight);
      _hist_y->fill(y,weight);
      _hist_x->fill(x,weight);

      /// @todo Do the event by event analysis here

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_hist_Q2); // normalize to unity
      normalize(_hist_y); // normalize to unity

    }

    //@}


    /// The histograms.
    Histo1DPtr _hist_Q2, _hist_y, _hist_x;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_DIS_Check);


}
