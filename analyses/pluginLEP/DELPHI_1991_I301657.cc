// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief DELPHI LEP1 charged multiplicity, code basically a copy of the ALEPH one
  /// @author Peter Richardson
  class DELPHI_1991_I301657 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1991_I301657);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      book(_histChTot, 2, 1, 1);
      book(_histAver , 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histChTot->fill(cfs.size());
      _histAver->fill(sqrtS(),cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histChTot, 200.0/sumOfWeights()); // bin width (2) and %age (100)
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histChTot;
    Profile1DPtr _histAver;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_1991_I301657);


}
