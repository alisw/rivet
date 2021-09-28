// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2011_S8971293 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2011_S8971293()
      : Analysis("ATLAS_2011_S8971293")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// Initialise and register projections
      declare(FastJets(FinalState(), FastJets::ANTIKT, 0.6), "AntiKtJets06");

      /// Book histograms
      {Histo1DPtr tmp; _h_deltaPhi.add(110., 160., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.add(160., 210., book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_deltaPhi.add(210., 260., book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_deltaPhi.add(260., 310., book(tmp, 1, 1, 4));}
      {Histo1DPtr tmp; _h_deltaPhi.add(310., 400., book(tmp, 1, 1, 5));}
      {Histo1DPtr tmp; _h_deltaPhi.add(400., 500., book(tmp, 1, 1, 6));}
      {Histo1DPtr tmp; _h_deltaPhi.add(500., 600., book(tmp, 1, 1, 7));}
      {Histo1DPtr tmp; _h_deltaPhi.add(600., 800., book(tmp, 1, 1, 8));}
      {Histo1DPtr tmp; _h_deltaPhi.add(800.,10000.,book(tmp, 1, 1, 9));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      Jets jets06;
      for (const Jet& jet : apply<FastJets>(event, "AntiKtJets06").jetsByPt(100.0*GeV)) {
        if (jet.absrap() < 2.8) {
          jets06.push_back(jet);
        }
      }
      if (jets06.size()>1){
        if (fabs(jets06[0].rapidity())<0.8 && fabs(jets06[1].rapidity())<0.8) {
          double observable = mapAngle0ToPi(jets06[0].phi()-jets06[1].phi()) / M_PI;
          _h_deltaPhi.fill(jets06[0].pT(), observable, weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr hist : _h_deltaPhi.histos()) {
        normalize(hist, 1/M_PI);
      }
    }

    //@}


  private:

    /// @name Histograms
    //@{
    BinnedHistogram _h_deltaPhi;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S8971293);

}
