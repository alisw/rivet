// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  class ATLAS_2011_I928289_W : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I928289_W()
      : Analysis("ATLAS_2011_I928289_W")
    {

    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ///Initialise and register projections here
      FinalState fs;

      Cut cut = (Cuts::pT >= 20*GeV);

      WFinder wfinder_el_bare(   fs, cut, PID::ELECTRON, 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.0, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      WFinder wfinder_el_dressed(fs, cut, PID::ELECTRON, 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.1, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      WFinder wfinder_mu_bare   (fs, cut, PID::MUON    , 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.0, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      WFinder wfinder_mu_dressed(fs, cut, PID::MUON    , 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.1, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);

      declare(wfinder_el_bare   , "WFinder_el_bare");
      declare(wfinder_el_dressed, "WFinder_el_dressed");
      declare(wfinder_mu_bare   , "WFinder_mu_bare");
      declare(wfinder_mu_dressed, "WFinder_mu_dressed");

      /// Book histograms here
      book(_h_Wminus_lepton_eta_el_bare       ,3, 1, 1);
      book(_h_Wminus_lepton_eta_el_dressed    ,3, 1, 2);
      book(_h_Wminus_lepton_eta_mu_bare       ,3, 1, 3);
      book(_h_Wminus_lepton_eta_mu_dressed    ,3, 1, 4);
      book(_h_Wplus_lepton_eta_el_bare        ,5, 1, 1);
      book(_h_Wplus_lepton_eta_el_dressed     ,5, 1, 2);
      book(_h_Wplus_lepton_eta_mu_bare        ,5, 1, 3);
      book(_h_Wplus_lepton_eta_mu_dressed     ,5, 1, 4);
      book(_h_W_asym_eta_el_bare              ,7, 1, 1);
      book(_h_W_asym_eta_el_dressed           ,7, 1, 2);
      book(_h_W_asym_eta_mu_bare              ,7, 1, 3);
      book(_h_W_asym_eta_mu_dressed           ,7, 1, 4);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const WFinder& wfinder_el_bare     = apply<WFinder>(event, "WFinder_el_bare");
      const WFinder& wfinder_el_dressed  = apply<WFinder>(event, "WFinder_el_dressed");
      const WFinder& wfinder_mu_bare     = apply<WFinder>(event, "WFinder_mu_bare");
      const WFinder& wfinder_mu_dressed  = apply<WFinder>(event, "WFinder_mu_dressed");

      fillPlots1D(wfinder_el_bare   , _h_Wplus_lepton_eta_el_bare   , _h_Wminus_lepton_eta_el_bare);
      fillPlots1D(wfinder_el_dressed, _h_Wplus_lepton_eta_el_dressed, _h_Wminus_lepton_eta_el_dressed);
      fillPlots1D(wfinder_mu_bare   , _h_Wplus_lepton_eta_mu_bare   , _h_Wminus_lepton_eta_mu_bare);
      fillPlots1D(wfinder_mu_dressed, _h_Wplus_lepton_eta_mu_dressed, _h_Wminus_lepton_eta_mu_dressed);
    }


    void fillPlots1D(const WFinder& wfinder, Histo1DPtr hist_plus, Histo1DPtr hist_minus) {
      if (wfinder.bosons().size() != 1) return;
      const Particle l = wfinder.constituentLeptons()[0];
      const FourMomentum miss = wfinder.constituentNeutrinos()[0];
      if (l.pT() > 20*GeV && miss.Et() > 25*GeV && wfinder.mT() > 40*GeV)
        (l.charge3() > 0 ? hist_plus : hist_minus)->fill(l.abseta());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
      divide(*_h_Wplus_lepton_eta_el_bare - *_h_Wminus_lepton_eta_el_bare,
             *_h_Wplus_lepton_eta_el_bare + *_h_Wminus_lepton_eta_el_bare,
             _h_W_asym_eta_el_bare);
      divide(*_h_Wplus_lepton_eta_el_dressed - *_h_Wminus_lepton_eta_el_dressed,
             *_h_Wplus_lepton_eta_el_dressed + *_h_Wminus_lepton_eta_el_dressed,
             _h_W_asym_eta_el_dressed);
      divide(*_h_Wplus_lepton_eta_mu_bare - *_h_Wminus_lepton_eta_mu_bare,
             *_h_Wplus_lepton_eta_mu_bare + *_h_Wminus_lepton_eta_mu_bare,
             _h_W_asym_eta_mu_bare);
      divide(*_h_Wplus_lepton_eta_mu_dressed - *_h_Wminus_lepton_eta_mu_dressed,
             *_h_Wplus_lepton_eta_mu_dressed + *_h_Wminus_lepton_eta_mu_dressed,
             _h_W_asym_eta_mu_dressed);

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());

      ///  Normalise, scale and otherwise manipulate histograms here
      const double sf = 0.5 * xs_pb / sumw; // 0.5 accounts for rapidity bin width
      scale(_h_Wminus_lepton_eta_el_bare   , sf);
      scale(_h_Wminus_lepton_eta_el_dressed, sf);
      scale(_h_Wminus_lepton_eta_mu_bare   , sf);
      scale(_h_Wminus_lepton_eta_mu_dressed, sf);
      scale(_h_Wplus_lepton_eta_el_bare    , sf);
      scale(_h_Wplus_lepton_eta_el_dressed , sf);
      scale(_h_Wplus_lepton_eta_mu_bare    , sf);
      scale(_h_Wplus_lepton_eta_mu_dressed , sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Wminus_lepton_eta_el_bare;
    Histo1DPtr _h_Wminus_lepton_eta_el_dressed;
    Histo1DPtr _h_Wminus_lepton_eta_mu_bare;
    Histo1DPtr _h_Wminus_lepton_eta_mu_dressed;
    Histo1DPtr _h_Wplus_lepton_eta_el_bare;
    Histo1DPtr _h_Wplus_lepton_eta_el_dressed;
    Histo1DPtr _h_Wplus_lepton_eta_mu_bare;
    Histo1DPtr _h_Wplus_lepton_eta_mu_dressed;
    Scatter2DPtr _h_W_asym_eta_el_bare;
    Scatter2DPtr _h_W_asym_eta_el_dressed;
    Scatter2DPtr _h_W_asym_eta_mu_bare;
    Scatter2DPtr _h_W_asym_eta_mu_dressed;

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I928289_W);

}
