// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF two-jet triply-differential cross-section
  class CDF_2001_S4517016 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2001_S4517016);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");

      {Histo1DPtr tmp; _h_ET.add(0.1, 0.7, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_ET.add(0.7, 1.4, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_ET.add(1.4, 2.1, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_ET.add(2.1, 3.0, book(tmp, 4, 1, 1));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 10*GeV, cmpMomByEt);
      if (jets.size() < 2) vetoEvent;
      FourMomentum jet1 = jets[0].momentum();
      FourMomentum jet2 = jets[1].momentum();
      double eta1 = jet1.abseta();
      double eta2 = jet2.abseta();
      double ET1 = jet1.Et();
      double ET2 = jet2.Et();
      if (!inRange(eta1, 0.1, 0.7) || ET1 < 40.0*GeV) vetoEvent;
      if (!inRange(eta2, 0.1, 3.0)) vetoEvent;
      _h_ET.fill(eta2, ET1);
      if (eta2<0.7 && ET2>40.0*GeV) _h_ET.fill(eta1, ET2);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double deta1 = 1.2;
      _h_ET.scale(crossSection()/nanobarn/sumOfWeights()/deta1 / 2.0, this);
    }

    //@}


  private:

    /// Histograms
    BinnedHistogram _h_ET;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2001_S4517016, CDF_2001_I538041);

}
