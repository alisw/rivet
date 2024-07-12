// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged multiplicity below Z pole based on ALEPH Z pole analysis
  // @author Peter Richardson
  class AMY_1990_I295160 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(AMY_1990_I295160);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      int offset = 0;
      if(isCompatibleWithSqrtS(50.0)) {
	offset = 1;
      }
      else if(isCompatibleWithSqrtS(52.0)) {
	offset = 2;
      }
      else if(isCompatibleWithSqrtS(55.0)) {
	offset = 3;
      }
      else if(isCompatibleWithSqrtS(56.0)) {
	offset = 4;
      }
      else if(isCompatibleWithSqrtS(57.0)) {
	offset = 5;
      }
      else if(isCompatibleWithSqrtS(60.0)) {
	offset = 6;
      }
      else if(isCompatibleWithSqrtS(60.8)) {
	offset = 7;
      }
      else if(isCompatibleWithSqrtS(61.4)) {
	offset = 8;
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book( _histChTotal , 1, 1, offset);
      book( _histTotal , 2, 1, 1);
      if(offset==5) {
	book( _histChAver, 1, 1, 9);
	book( _histAver, 2, 2, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histChTotal->fill(cfs.size());
      _histTotal->fill(sqrtS()/GeV,cfs.size());
      if(_histAver) {
	_histChAver->fill(cfs.size());
	_histAver->fill(sqrtS()/GeV,cfs.size());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      scale(_histChTotal, 200.0/sumOfWeights()); // bin width (2) and %age (100)
      if(_histAver) 
	scale(_histChAver, 200.0/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histChTotal;
    Histo1DPtr _histChAver;
    Profile1DPtr _histTotal;
    Profile1DPtr _histAver;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(AMY_1990_I295160);


}
