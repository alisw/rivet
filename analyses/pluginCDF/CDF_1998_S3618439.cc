// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF diff cross-section in events with large missing energy
  class CDF_1998_S3618439 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1998_S3618439);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");

      book(_h_sumET_20 ,1, 1, 1);
      book(_h_sumET_100 ,1, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 20*GeV, cmpMomByEt);
      double sumET_20(0.0), sumET_100(0.0);
      for (const Jet& jet : jets) {
        double ET = jet.Et()/GeV;
        sumET_20 += ET;
        if (ET > 100.0) sumET_100 += ET;
      }
      if (sumET_20 > 320.0) _h_sumET_20->fill(sumET_20);
      if (sumET_100 > 320.0) _h_sumET_100->fill(sumET_100);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sumET_20, crossSection()/picobarn/sumOfWeights());
      scale(_h_sumET_100, crossSection()/picobarn/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_sumET_20, _h_sumET_100;
    //@}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_1998_S3618439, CDF_1998_I448075);

}
