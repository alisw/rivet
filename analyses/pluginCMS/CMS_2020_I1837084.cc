// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Measurements of differential Z-boson -> ll and vv production cross-sections in 13 TeV pp collisions
  class CMS_2020_I1837084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2020_I1837084);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      ZFinder zmmFind(FinalState(), Cuts::pT > 0*GeV, PID::MUON, 76.1876*GeV, 106.1876*GeV, 0.1,
                      ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES );
      declare(zmmFind, "ZmmFind");

      // Book histograms
      book(_h_Z_pt,      12, 1, 1);
      book(_h_Z_pt_norm, 13, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& zmms = apply<ZFinder>(event, "ZmmFind").bosons();

      if (zmms.size() == 1 && zmms[0].pT() > 200*GeV) {
        _h_Z_pt     ->fill(min(zmms[0].pT()/GeV, 1499.999));
        _h_Z_pt_norm->fill(min(zmms[0].pT()/GeV, 1499.999));
      }

    }


    /// @todo Replace with barchart()
    void normalizeToSum(Histo1DPtr hist) {
      double sum = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        sum += hist->bin(i).height();
        double width = hist->bin(i).width();
        hist->bin(i).scaleW(width != 0 ? width : 1.);
      }
      if (hist->integral() > 0) scale(hist, 1./hist->integral());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/femtobarn/sumOfWeights() : 1.0;

      scale(_h_Z_pt, norm);

      normalizeToSum(_h_Z_pt_norm);

    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Z_pt, _h_Z_pt_norm;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2020_I1837084);

}
