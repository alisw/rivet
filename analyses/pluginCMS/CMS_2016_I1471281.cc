// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Measurement of the W and Z bosons pT produced in pp collisions at sqrt(s)=8 TeV
  class CMS_2016_I1471281 : public Analysis {
    public:

      /// Constructor
      CMS_2016_I1471281(std::string name="CMS_2016_I1471281")
        : Analysis(name)
        {
          _mode = 0; // init 
        }


      /// @name Analysis methods
      //@{

      /// Book histograms and initialise projections before the run
      void init() {
        
        // Get options from the new option system
        // default to both.
        if ( getOption("VMODE") == "BOTH" ) _mode = 0;
        if ( getOption("VMODE") == "W" )    _mode = 1;
        if ( getOption("VMODE") == "Z" )    _mode = 2;

        // Set up projections
        FinalState fs;

        Cut cut_mu = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;

        // Dressed Ws ...
        WFinder wmunu_Finder(fs, cut_mu, PID::MUON, 0*GeV, YODA::MAXDOUBLE, 0*GeV, 0, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
        declare(wmunu_Finder, "Wmunu_Finder");

        // Dressed Zs ... 
        ZFinder zmumu_Finder(fs, cut_mu, PID::MUON, 60*GeV, 120*GeV, 0, ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
        declare(zmumu_Finder, "Zmumu_Finder");

        // Histograms
        if (_mode == 0 || _mode == 1) {
          book(_hist_WtoMuNuPt, 8, 1, 1);
        }
        if (_mode == 0 || _mode == 2) {
          book(_hist_ZtoMuMuPt, 9, 1, 1);
        }
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        if (_mode == 0 || _mode == 1) {
          // Get the W bosons - muon decay channel
          const WFinder& wmunu_Finder = apply<WFinder>(event, "Wmunu_Finder");
          if (!wmunu_Finder.bosons().empty()) {
            const FourMomentum pWmunu = wmunu_Finder.bosons()[0].momentum();
            _hist_WtoMuNuPt->fill(pWmunu.pT()/GeV);
          }
        }

        if (_mode == 0 || _mode == 2) {
          // Get the Z bosons - muon decay channel
          const ZFinder& zmumu_Finder = apply<ZFinder>(event, "Zmumu_Finder");
          if (!zmumu_Finder.bosons().empty()) {
            const FourMomentum pZmumu = zmumu_Finder.bosons()[0].momentum();
            _hist_ZtoMuMuPt->fill(pZmumu.pT()/GeV);
          }
        }

      }


      /// Normalise histograms etc., after the run
      void finalize() {

        MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
        MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
        MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());

        if (_mode == 0 || _mode == 1) {
          normalize(_hist_WtoMuNuPt);
        }

        if (_mode == 0 || _mode == 2) {
          normalize(_hist_ZtoMuMuPt);
        }
              

      }

      //@}

    protected:

      // Data members
      size_t _mode;


    private:


      /// @name Histograms

      Histo1DPtr _hist_WtoMuNuPt;
      Histo1DPtr _hist_ZtoMuMuPt;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1471281);

}
