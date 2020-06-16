// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class ATLAS_2014_I1300647 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1300647()
      : Analysis("ATLAS_2014_I1300647")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
   
      FinalState fs;
      
      ZFinder zfinder_dressed_el(fs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1);
      declare(zfinder_dressed_el, "ZFinder_dressed_el");
      
      ZFinder zfinder_bare_el(fs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.0);
      declare(zfinder_bare_el,	"ZFinder_bare_el");
      
      ZFinder zfinder_dressed_mu(fs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::MUON,     66.0*GeV, 116.0*GeV, 0.1);
      declare(zfinder_dressed_mu, "ZFinder_dressed_mu");
      
      ZFinder zfinder_bare_mu(fs, Cuts::abseta<2.4 && Cuts::pT>20.0*GeV, PID::MUON,     66.0*GeV, 116.0*GeV, 0.0);
      declare(zfinder_bare_mu,	"ZFinder_bare_mu");
      
      // Book histograms
      book(_hist_zpt_el_dressed ,1, 1, 1);  // electron "dressed"
      book(_hist_zpt_mu_dressed ,1, 1, 2);  // muon "dressed"
      book(_hist_zpt_el_bare    ,1, 2, 1);  // electron "bare"
      book(_hist_zpt_mu_bare    ,1, 2, 2);  // muon "bare"	 

      //double-differential plots
      {Histo1DPtr tmp; _h_zpt_el_mu_dressed.add(0.0, 1.0, book(tmp, 3, 1, 2));}
      {Histo1DPtr tmp; _h_zpt_el_mu_dressed.add(1.0, 2.0, book(tmp, 3, 1, 4));}
      {Histo1DPtr tmp; _h_zpt_el_mu_dressed.add(2.0, 2.4, book(tmp, 3, 1, 6));}
 
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder_dressed_el = apply<ZFinder>(event, "ZFinder_dressed_el");
      const ZFinder& zfinder_bare_el    = apply<ZFinder>(event, "ZFinder_bare_el");
      const ZFinder& zfinder_dressed_mu = apply<ZFinder>(event, "ZFinder_dressed_mu");
      const ZFinder& zfinder_bare_mu    = apply<ZFinder>(event, "ZFinder_bare_mu");	
      
      FillPlots1d(zfinder_dressed_el, _hist_zpt_el_dressed);

      FillPlots1d(zfinder_bare_el,    _hist_zpt_el_bare);
          
      FillPlots1d(zfinder_dressed_mu, _hist_zpt_mu_dressed);
          
      FillPlots1d(zfinder_bare_mu,    _hist_zpt_mu_bare);  
      
      FillPlots3d(zfinder_dressed_el, _h_zpt_el_mu_dressed);      
      FillPlots3d(zfinder_dressed_mu, _h_zpt_el_mu_dressed);    

    }

    void FillPlots1d(const ZFinder& zfinder, Histo1DPtr hist) {
      if(zfinder.bosons().size() != 1) return;
      const FourMomentum pZ = zfinder.bosons()[0].momentum();
      hist->fill(pZ.pT()/GeV);
      return; 
    } 
    
    void FillPlots3d(const ZFinder& zfinder, BinnedHistogram& binnedHist) {
      if(zfinder.bosons().size() != 1) return;
      const FourMomentum pZ = zfinder.bosons()[0].momentum();
      binnedHist.fill(pZ.rapidity(), pZ.pT()/GeV);   
      return; 
    }  

    /// Normalise histograms etc., after the run
    void finalize() {
    
      normalize(_hist_zpt_el_dressed);
      normalize(_hist_zpt_el_bare);
      
      normalize(_hist_zpt_mu_dressed);  
      normalize(_hist_zpt_mu_bare); 

      for (Histo1DPtr hist : _h_zpt_el_mu_dressed.histos()) { normalize(hist); }

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    BinnedHistogram _h_zpt_el_mu_dressed;
  
 
    Histo1DPtr _hist_zpt_el_dressed;
    Histo1DPtr _hist_zpt_el_bare;     
    Histo1DPtr _hist_zpt_mu_dressed;
    Histo1DPtr _hist_zpt_mu_bare;    

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1300647);
}
