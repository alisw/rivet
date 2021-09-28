// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#include "Rivet/Projections/TriggerCDFRun0Run1.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// @brief E735 charged multiplicity in NSD-triggered events
  class E735_1998_S3905616 : public Analysis {
  public:

    /// Constructor
    E735_1998_S3905616() : Analysis("E735_1998_S3905616") {}


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState(), "FS");

      // Histo
      book(_hist_multiplicity ,1, 1, 1);
      book(_sumWTrig, "TMP/sumWtrig");
    
    }


    void analyze(const Event& event) {
      const bool trigger = apply<TriggerUA5>(event, "Trigger").nsdDecision();
      if (!trigger) vetoEvent;
      _sumWTrig->fill();

      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      _hist_multiplicity->fill(numParticles);
    }


    void finalize() {
      scale(_hist_multiplicity, 1 / *_sumWTrig);
    }

    //@}


  private:

    /// @name Weight counter
    //@{
    CounterPtr _sumWTrig;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _hist_multiplicity;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(E735_1998_S3905616);

}
