// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /// Electroweak Wjj production at 8 TeV
  class ATLAS_2013_I1217863 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1217863);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 2;
      _doZ  = true;
      _doW  = true;
      if ( getOption("LMODE") == "EL" ) { _mode = 2;}
      if ( getOption("LMODE") == "MU" ) _mode = 3;
      if ( getOption("LMODE") == "ZEL" ) {
        _mode = 2;
        _doW  = false;
      }
      if ( getOption("LMODE") == "ZMU" ) {
        _mode = 3;
        _doW  = false;
      }
      if ( getOption("LMODE") == "WEL" ) {
        _mode = 2;
        _doZ  = false;
      }
      if ( getOption("LMODE") == "WMU" ) {
        _mode = 3;
        _doZ  = false;
      }

      FinalState fs;
      declare(fs, "FS");

      Cut cuts = Cuts::abseta < 2.47 && Cuts::pT > 25*GeV;

      // Z finder
      if (_doZ) {
        ZFinder zf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 40.0*GeV, 1000.0*GeV, 0.1, 
                   ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
        declare(zf, "ZF");
      }

      if (_doW) {
        // W finder for electrons and muons
        WFinder wf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 0.0*GeV, 1000.0*GeV, 35.0*GeV, 0.1,
                   WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
        declare(wf, "WF");
      }

      // leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.37 && Cuts::pT > 15*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // jets
      VetoedFinalState jet_fs(fs);
      if (_doZ) { jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF")); }
      if (_doW) { jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF")); }
      jet_fs.addVetoOnThisFinalState(getProjection<LeadingParticlesFinalState>("LeadingPhoton"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      declare(jets, "Jets");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "isolatedFS");


      // Book histograms
      if (_doZ) {
        book(_hist_EgammaT_inclZ   ,11, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
        book(_hist_EgammaT_exclZ   ,12, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
        book(_hist_Njet_EgammaT15Z ,17, 1, _mode); // dSigma / dNjet for E^gamma_T >= 15
        book(_hist_Njet_EgammaT60Z ,18, 1, _mode); // dSigma / dNjet for E^gamma_T >= 60
        book(_hist_mZgamma         ,20, 1, _mode); // dSigma / dm^{Zgamma}
      }
      if (_doW){
        book(_hist_EgammaT_inclW   , 7, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
        book(_hist_EgammaT_exclW   , 8, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
        book(_hist_Njet_EgammaT15W ,15, 1, _mode); // dSigma / dNjet for E^gamma_T >= 15
        book(_hist_Njet_EgammaT60W ,16, 1, _mode); // dSigma / dNjet for E^gamma_T >= 60
        book(_hist_mWgammaT        ,19, 1, _mode); // dSigma / dm^{Zgamma}
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // retrieve leading photon
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1)  vetoEvent;
      const Particle& leadingPhoton = photons[0];
      if (leadingPhoton.Et() < 15.0*GeV) vetoEvent;
      if (leadingPhoton.abseta() > 2.37) vetoEvent;

      // check photon isolation
      double coneEnergy(0.0);
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      for (const Particle& p : fs) {
        if ( deltaR(leadingPhoton, p) < 0.4 )  coneEnergy += p.E();
      }
      if (coneEnergy / leadingPhoton.E() >= 0.5 )  vetoEvent;

      if (_doW) {
	// retrieve W boson candidate
	const WFinder& wf = apply<WFinder>(event, "WF");
	if ( wf.bosons().size() == 1 ) { 
	  
	  // retrieve constituent neutrino
	  const Particle& neutrino = wf.constituentNeutrino();
	  if ( (neutrino.pT() > 35.0*GeV) ) {
	    
	    // retrieve constituent lepton
	    const Particle& lepton = wf.constituentLepton();
	    if ( lepton.pT() > 25.0*GeV && lepton.abseta() < 2.47 ) {
	      
	      // check photon-lepton overlap
	      if ( deltaR(leadingPhoton, lepton) > 0.7 ) {
		
		// count jets
		const FastJets& jetfs = apply<FastJets>(event, "Jets");
		Jets jets = jetfs.jets(cmpMomByEt);
		int goodJets = 0;
		for (const Jet& j : jets) {
		  if ( !(j.Et() > 30.0*GeV) )  break;
		  if ( (j.abseta() < 4.4) &&				\
		       (deltaR(leadingPhoton, j) > 0.3) &&		\
		       (deltaR(lepton,        j) > 0.3) )  ++goodJets;
		}
		
		double Njets = double(goodJets) + 0.5;
		double photonEt = leadingPhoton.Et()*GeV;
		
		const FourMomentum& lep_gamma = lepton.momentum() + leadingPhoton.momentum();
		double term1 = sqrt(lep_gamma.mass2() + lep_gamma.pT2()) + neutrino.Et();
		double term2 = (lep_gamma + neutrino.momentum()).pT2();
		double mWgammaT = sqrt(term1 * term1 - term2) * GeV;
		
		_hist_EgammaT_inclW->fill(photonEt);
		
		_hist_Njet_EgammaT15W->fill(Njets);
		
		if ( !goodJets )  _hist_EgammaT_exclW->fill(photonEt);
		
		if (photonEt > 40.0*GeV) {
		  _hist_mWgammaT->fill(mWgammaT);
		  if (photonEt > 60.0*GeV)  _hist_Njet_EgammaT60W->fill(Njets);
		}
	      }
	    }
	  }
	}
      }

      if (_doZ ){

	// retrieve Z boson candidate
	const ZFinder& zf = apply<ZFinder>(event, "ZF");
	if ( zf.bosons().size() == 1 ) {
	  const Particle& Zboson  = zf.boson();
	  if ( (Zboson.mass() > 40.0*GeV) ) {
	    
	    // check charge of constituent leptons
	    const Particles& leptons = zf.constituents();
	    if (leptons.size() == 2 && leptons[0].charge() * leptons[1].charge() < 0.) {
	      
	      bool lpass = true;
	      // check photon-lepton overlap
	      for (const Particle& p : leptons) {
		if ( !(p.pT() > 25.0*GeV && p.abseta() < 2.47 && deltaR(leadingPhoton, p) > 0.7) )  lpass = false;
	      }
	      if ( lpass ) {
		
		// count jets
		const FastJets& jetfs = apply<FastJets>(event, "Jets");
		Jets jets = jetfs.jets(cmpMomByEt);
		int goodJets = 0;
		for (const Jet& j : jets) {
		  if ( !(j.Et() > 30.0*GeV) )  break;
		  if ( (j.abseta() < 4.4) &&		    \
		       (deltaR(leadingPhoton, j) > 0.3) &&  \
		       (deltaR(leptons[0],    j) > 0.3) &&		\
		       (deltaR(leptons[1],    j) > 0.3) )  ++goodJets;
		}
		
		double Njets = double(goodJets) + 0.5;
		double photonEt = leadingPhoton.Et()*GeV;
		double mZgamma = (Zboson.momentum() + leadingPhoton.momentum()).mass() * GeV;
		
		_hist_EgammaT_inclZ->fill(photonEt);
		
		_hist_Njet_EgammaT15Z->fill(Njets);
		
		if ( !goodJets )   _hist_EgammaT_exclZ->fill(photonEt);
		
		if (photonEt >= 40.0*GeV) {
		  _hist_mZgamma->fill(mZgamma);
		  if (photonEt >= 60.0*GeV)  _hist_Njet_EgammaT60Z->fill(Njets);
		}
	      }
	    }
	  }
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_fb = crossSection()/femtobarn;
      const double sumw = sumOfWeights();
      const double sf = xs_fb / sumw;

      if (_doZ) {
	scale(_hist_EgammaT_exclZ, sf);
	scale(_hist_EgammaT_inclZ, sf);
	normalize(_hist_Njet_EgammaT15Z);
	normalize(_hist_Njet_EgammaT60Z);
	normalize(_hist_mZgamma);
      }

      if (_doW) {
	scale(_hist_EgammaT_exclW, sf);
	scale(_hist_EgammaT_inclW, sf);
	normalize(_hist_Njet_EgammaT15W);
	normalize(_hist_Njet_EgammaT60W);
	normalize(_hist_mWgammaT);
      }

    }

    //@}

  protected:

    size_t _mode;
    bool _doW;
    bool _doZ;

  private:

    /// @name Histograms
    //@{

    Histo1DPtr _hist_EgammaT_inclZ;
    Histo1DPtr _hist_EgammaT_exclZ;
    Histo1DPtr _hist_Njet_EgammaT15Z;
    Histo1DPtr _hist_Njet_EgammaT60Z;
    Histo1DPtr _hist_mZgamma;

    Histo1DPtr _hist_EgammaT_inclW;
    Histo1DPtr _hist_EgammaT_exclW;
    Histo1DPtr _hist_Njet_EgammaT15W;
    Histo1DPtr _hist_Njet_EgammaT60W;
    Histo1DPtr _hist_mWgammaT;

    //@}

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1217863);

}

