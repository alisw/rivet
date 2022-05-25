// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1989_I277658 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1989_I277658);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      int offset = 0;
      if(isCompatibleWithSqrtS(14.0)) {
	offset = 1;
      }
      else if(isCompatibleWithSqrtS(22.0)) {
	offset = 2;
      }
      else if(isCompatibleWithSqrtS(34.8)) {
	offset = 3;
      }
      else if(isCompatibleWithSqrtS(43.6)) {
	offset = 4;
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book(_histCh, 5, 1, offset); 
      book(_histTotal, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histCh->fill(cfs.size());
      _histTotal->fill(sqrtS(),cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histCh, 2.0/sumOfWeights()); // bin width (2)
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histCh;
    Profile1DPtr _histTotal;
    //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1989_I277658);


}
