// -*- C++ -*-
#include <cmath>

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"


namespace Rivet {

  using namespace Cuts;


  class ATLAS_2011_I928289_W : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I928289_W()
      : Analysis("ATLAS_2011_I928289_W")
    {
      setNeedsCrossSection(true);
    }


  public:

    /// @name Analysis methods
    //@{
    
    /// Book histograms and initialise projections before the run
    void init() {

      ///Initialise and register projections here
      FinalState fs;
      
      Cut cut = pT >= 20*GeV;

      WFinder wfinder_el_bare(   fs, cut, PID::ELECTRON, 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.0, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      WFinder wfinder_el_dressed(fs, cut, PID::ELECTRON, 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      WFinder wfinder_mu_bare   (fs, cut, PID::MUON    , 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.0, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      WFinder wfinder_mu_dressed(fs, cut, PID::MUON    , 40.0*GeV, 7000.0*GeV, 25.0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);

      addProjection(wfinder_el_bare   , "WFinder_el_bare");
      addProjection(wfinder_el_dressed, "WFinder_el_dressed");
      addProjection(wfinder_mu_bare   , "WFinder_mu_bare");
      addProjection(wfinder_mu_dressed, "WFinder_mu_dressed");

      /// Book histograms here
      _h_Wminus_lepton_eta_el_bare       = bookHisto1D(3, 1, 1);
      _h_Wminus_lepton_eta_el_dressed    = bookHisto1D(3, 1, 2);
      _h_Wminus_lepton_eta_mu_bare       = bookHisto1D(3, 1, 3);
      _h_Wminus_lepton_eta_mu_dressed    = bookHisto1D(3, 1, 4);
      _h_Wplus_lepton_eta_el_bare        = bookHisto1D(5, 1, 1);
      _h_Wplus_lepton_eta_el_dressed     = bookHisto1D(5, 1, 2);
      _h_Wplus_lepton_eta_mu_bare        = bookHisto1D(5, 1, 3);
      _h_Wplus_lepton_eta_mu_dressed     = bookHisto1D(5, 1, 4);
      _h_W_asym_eta_el_bare              = bookScatter2D(7, 1, 1);
      _h_W_asym_eta_el_dressed           = bookScatter2D(7, 1, 2);
      _h_W_asym_eta_mu_bare              = bookScatter2D(7, 1, 3);
      _h_W_asym_eta_mu_dressed           = bookScatter2D(7, 1, 4);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      ///Do the event by event analysis here
      const WFinder& wfinder_el_bare     = applyProjection<WFinder>(event, "WFinder_el_bare");
      const WFinder& wfinder_el_dressed  = applyProjection<WFinder>(event, "WFinder_el_dressed");
      const WFinder& wfinder_mu_bare     = applyProjection<WFinder>(event, "WFinder_mu_bare");
      const WFinder& wfinder_mu_dressed  = applyProjection<WFinder>(event, "WFinder_mu_dressed");

      FillPlots1d(wfinder_el_bare   , _h_Wplus_lepton_eta_el_bare   , _h_Wminus_lepton_eta_el_bare   , weight);
      FillPlots1d(wfinder_el_dressed, _h_Wplus_lepton_eta_el_dressed, _h_Wminus_lepton_eta_el_dressed, weight);
      FillPlots1d(wfinder_mu_bare   , _h_Wplus_lepton_eta_mu_bare   , _h_Wminus_lepton_eta_mu_bare   , weight);
      FillPlots1d(wfinder_mu_dressed, _h_Wplus_lepton_eta_mu_dressed, _h_Wminus_lepton_eta_mu_dressed, weight);


    }

    void FillPlots1d(const WFinder& wfinder, Histo1DPtr hist_plus, Histo1DPtr hist_minus, double weight) {

      if (wfinder.bosons().size() != 1) return;

      Particle l = wfinder.constituentLeptons()[0];
      const FourMomentum& miss = wfinder.constituentNeutrinos()[0].momentum();

      if(l.momentum().pT() > 20*GeV && miss.Et() > 25*GeV && wfinder.mT() > 40*GeV) {
        int lepCharge = l.charge(); 
        if      (lepCharge > 0) hist_plus ->fill( fabs(l.eta()), weight);
        else if (lepCharge < 0) hist_minus->fill( fabs(l.eta()), weight);
      }

      return; 
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
      MSG_INFO( "Cross-Section/pb     : " << xs_pb       );
      MSG_INFO( "Sum of weights       : " << sumw        );
      MSG_INFO( "nEvents              : " << numEvents() );

      const double sf(0.5 * xs_pb / sumw); // 0.5 accounts for rapidity bin width

      ///  Normalise, scale and otherwise manipulate histograms here
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

    // Data members like post-cuts event weight counters go here


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
