// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @brief inclusive W/Z cross sections at 7 TeV
  class ATLAS_2016_I1502620 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1502620);
    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // run everything
      _mode = 0;
      _runZ = true;
      _runW = true;
      if ( getOption("LMODE") == "EL" || 
	   getOption("LMODE") == "ZEL" ||
	   getOption("LMODE") == "WEL" ) 
	_mode = 1;
      if ( getOption("LMODE") == "MU" || 
	   getOption("LMODE") == "ZMU" ||
	   getOption("LMODE") == "WMU" ) 
	_mode = 2;
      if ( getOption("LMODE") == "Z" || 
	   getOption("LMODE") == "ZEL" || 
	   getOption("LMODE") == "ZMU" ) 
	_runW = false;
      if ( getOption("LMODE") == "W" || 
	   getOption("LMODE") == "WEL" || 
	   getOption("LMODE") == "WMU" ) 
	_runZ = false;




      ///Initialise and register projections here
      const FinalState fs;

      Cut Wcuts = Cuts::pT >= 25*GeV; // minimum lepton pT
      Cut Zcuts = Cuts::pT >= 20.0*GeV;

      WFinder wfinder_edressed(fs, Wcuts, PID::ELECTRON, 40*GeV, 13*TeV, 25*GeV, 0.1, 
				 WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wfinder_edressed, "WFinder_edressed");

      ZFinder zfindere(fs, Zcuts, PID::ELECTRON, 46.0*GeV, 150*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfindere, "ZFindere");

      WFinder wfinder_mdressed(fs, Wcuts, PID::MUON, 40*GeV, 13*TeV, 25*GeV, 0.1, 
				 WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wfinder_mdressed, "WFinder_mdressed");

      ZFinder zfinderm(fs, Zcuts, PID::MUON, 46.0*GeV, 150*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinderm, "ZFinderm");


      /// Book histograms here      
      if (_runW) {
	_h_Wp_eta = bookHisto1D(   9, 1, 1);
	_h_Wm_eta = bookHisto1D(  10, 1, 1);
	_h_W_asym = bookScatter2D(35, 1, 1);
      }

      if (_runZ) {
	_h_Zcenlow_y_dressed   = bookHisto1D(11, 1, 1);
	_h_Zcenpeak_y_dressed  = bookHisto1D(12, 1, 1);
	_h_Zcenhigh_y_dressed  = bookHisto1D(13, 1, 1);
	_h_Zfwdpeak_y_dressed  = bookHisto1D(14, 1, 1);
	_h_Zfwdhigh_y_dressed  = bookHisto1D(15, 1, 1);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      // W stuff 
      const WFinder& wfindere = apply<WFinder>(event, "WFinder_edressed");	     
      const WFinder& wfinderm = apply<WFinder>(event, "WFinder_mdressed");	     
      
      if (wfindere.bosons().size()+wfinderm.bosons().size() == 1 && _runW) {

        const double weight = event.weight();

	Particle lep;
	if (_mode !=2 && wfindere.bosons().size() == 1 ) {
	  lep = wfindere.constituentLeptons()[0];
	}
	else if (_mode !=1 && wfinderm.bosons().size() == 1 ) {
	  lep = wfinderm.constituentLeptons()[0];
	}
	if (lep.charge3() > 0)  _h_Wp_eta->fill(lep.abseta(), weight);
	else                    _h_Wm_eta->fill(lep.abseta(), weight);
	
      }

      // now the Z stuff. 
      const ZFinder& zfindere = apply<ZFinder>(event, "ZFindere");
      const ZFinder& zfinderm = apply<ZFinder>(event, "ZFinderm");
      

      // must be one and only one candidate.
      if (zfindere.bosons().size()+zfinderm.bosons().size() == 1 && _runZ) {

	const double weight = event.weight();
	
	Particle Zboson;
	ParticleVector leptons;
	
	// candidate is e+e-
	if (_mode != 2 && zfindere.bosons().size() == 1 ) {
	  
	  Zboson = zfindere.boson();
	  leptons = zfindere.constituents();
	}  

	// candidate is mu+mu-
        else if (_mode !=1 && zfinderm.bosons().size() == 1 ) {
	  
	  Zboson = zfinderm.boson();
	  leptons = zfinderm.constituents();
	  
	}
	const double zrap  = Zboson.absrap();
	const double zmass = Zboson.mass();
	const double eta1 = leptons[0].abseta();
	const double eta2 = leptons[1].abseta();
		
	// separation into central/forward and three mass bins
	if (eta1 < 2.5 && eta2 < 2.5) {
	  if (zmass < 66.0*GeV)        _h_Zcenlow_y_dressed->fill(zrap, weight);
	  else if (zmass < 116.0*GeV)  _h_Zcenpeak_y_dressed->fill(zrap, weight);
	  else                         _h_Zcenhigh_y_dressed->fill(zrap, weight);
	} 
	else if ((eta1 < 2.5 && 2.5 < eta2 && eta2 < 4.9) || (eta2 < 2.5 && 2.5 < eta1 && eta1 < 4.9)) {
	  if (zmass < 66.0*GeV)   vetoEvent;
	  if (zmass < 116.0*GeV)  _h_Zfwdpeak_y_dressed->fill(zrap, weight);
	  else                    _h_Zfwdhigh_y_dressed->fill(zrap, weight);
	}
      }
      
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
      //divide(*_h_Wp_eta - *_h_Wm_eta, *_h_Wp_eta + *_h_Wm_eta, _h_W_asym);
      if (_runW) {
	for (size_t i = 0; i < _h_Wp_eta->numBins(); ++i) {
	  YODA::HistoBin1D& bp = _h_Wp_eta->bin(i);
	  YODA::HistoBin1D& bm = _h_Wm_eta->bin(i);
	  const double sum  = bp.height() + bm.height();
	  //const double xerr = 0.5 * bp.xWidth();
	  double val = 0., yerr = 0.;

	  if (sum) {
	    const double pos2  = bp.height() * bp.height();
	    const double min2  = bm.height() * bm.height();
	    const double errp2 = bp.heightErr() * bp.heightErr();
	    const double errm2 = bm.heightErr() * bm.heightErr();
	    val = (bp.height() - bm.height()) / sum;
	    yerr = 2. * sqrt(errm2 * pos2 + errp2 * min2) / (sum * sum);
	  }
	  _h_W_asym->addPoint(bp.midpoint(), val, 0.5*bp.xWidth(), yerr);
	}
      }

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());
      MSG_DEBUG( "Cross-section/pb     : " << xs_pb       );
      MSG_DEBUG( "Sum of weights       : " << sumw        );
      MSG_DEBUG( "nEvents              : " << numEvents() );

      ///  Normalise, scale and otherwise manipulate histograms here
      double lfac = 1.0;
      // If we have been running on both electrons and muons, need to divide by two to
      // get the xsec for one flavour
      if (_mode == 0) lfac = 0.5;
      const double sf = lfac * 0.5 * xs_pb / sumw; // 0.5 accounts for rapidity bin width

      if (_runW){
	scale(_h_Wp_eta, sf);
	scale(_h_Wm_eta, sf);
      }

      if (_runZ){
	scale(_h_Zcenlow_y_dressed, sf);
	scale(_h_Zcenpeak_y_dressed, sf);
	scale(_h_Zcenhigh_y_dressed, sf);
	scale(_h_Zfwdpeak_y_dressed, sf);
	scale(_h_Zfwdhigh_y_dressed, sf);
      }
    }

    //@}


  protected:
    size_t _mode;
    bool _runZ, _runW;

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Wp_eta, _h_Wm_eta;
    Scatter2DPtr _h_W_asym;

    Histo1DPtr _h_Zcenlow_y_dressed;
    Histo1DPtr _h_Zcenpeak_y_dressed;
    Histo1DPtr _h_Zcenhigh_y_dressed;
    Histo1DPtr _h_Zfwdpeak_y_dressed;
    Histo1DPtr _h_Zfwdhigh_y_dressed;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620);


}

// END END END




