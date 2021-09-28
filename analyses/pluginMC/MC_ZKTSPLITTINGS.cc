// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  


  /// @brief MC validation analysis for Z + jets events
  class MC_ZKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_ZKTSPLITTINGS()
      : MC_JetSplittings("MC_ZKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
		  _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      FinalState fs;
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zfinder(fs, cut, _lepton, 65*GeV, 115*GeV, _dR, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::KT, 0.6);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;

      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    //@}

  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZKTSPLITTINGS);
}
