// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief EW Zjj using early Run-2 data
  class ATLAS_2017_I1627873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2017_I1627873);


    /// Book histograms and initialise projections before the run
    void init() {
    
      _mode = 0;
      if ( getOption("TYPE") == "EW_ONLY" ) _mode = 1;
      
      FinalState fs(Cuts::abseta < 5.0);

      FinalState photon_fs(Cuts::abspid == PID::PHOTON);

      PromptFinalState electron_fs(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muon_fs(Cuts::abspid == PID::MUON);

      DressedLeptons dressed_electrons(photon_fs, electron_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_electrons, "DressedElectrons");

      DressedLeptons dressed_muons(photon_fs, muon_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_muons, "DressedMuons");
			
      VetoedFinalState remfs(fs);
      remfs.addVetoOnThisFinalState(dressed_electrons);
      remfs.addVetoOnThisFinalState(dressed_muons);

      FastJets jets(remfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets, "Jets");
      
      if (_mode)  book(_h["zjj-ew"], 3, 1, 1);
      else        book(_h["zjj"], 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.4);
	    vector<DressedLepton> electrons = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons();
	    vector<DressedLepton> muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons();
	  
	   	// Overlap Removal
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);
      idiscardIfAnyDeltaRLess(muons,     jets, 0.4);

      Particle lep1, lep2;
      if (electrons.size() == 2 && muons.empty()) {
				lep1 = electrons[0]; lep2 = electrons[1];
				if ( lep1.charge3() == lep2.charge3() )  vetoEvent;
			} 
      else if (electrons.empty() && muons.size() == 2) {
				lep1 = muons[0]; lep2 = muons[1];
				if (lep1.charge3() == lep2.charge3())  vetoEvent;
			} 
      else  vetoEvent;
	   
      if (jets.size() < 2) vetoEvent;
      
      const FourMomentum dilepton = lep1.mom()+lep2.mom();
      if ( !inRange(dilepton.mass(), 81.0*GeV, 101.0*GeV) ) vetoEvent; 

      const double jet1pt = jets[0].pT();
      const double jet2pt = jets[1].pT();
      const double  mjj = (jets[0].mom() + jets[1].mom()).mass();
      const double  zpt = (lep1.mom() + lep2.mom()).pT();
      
      size_t ngapjets = 0;
      Jet thirdjet;
      for (size_t i = 2; i < jets.size(); ++i) { 
        const Jet j = jets[i];
        if (_isBetween(j, jets[0], jets[1])) {
          if (!ngapjets)  thirdjet = j;
          ++ngapjets;
        }
      } 
      
      const double ptbal_vec = (jets[0].mom() + jets[1].mom() + lep1.mom() + lep2.mom()).pT();
      const double ptbal_sc = jets[0].pT() + jets[1].pT() + lep1.pT() + lep2.pT();
      const double ptbalance2 = ptbal_vec / ptbal_sc;

      const double ptbal3_vec = (jets[0].mom() + jets[1].mom() + thirdjet.mom() + lep1.mom() + lep2.mom()).pT();
      const double ptbal3_sc = jets[0].pT() + jets[1].pT() + thirdjet.pT() + lep1.pT() + lep2.pT();
      const double ptbalance3 = ptbal3_vec / ptbal3_sc;



      //categories: baseline, high-PT, EW-enriched, QCD-enriched, high-mass, EW-enriched and high-mass
      if(!(jet1pt > 55*GeV && jet2pt > 45*GeV))  vetoEvent;

      if (_mode) { 
        if (zpt > 20.0*GeV && !ngapjets && ptbalance2 < 0.15 && mjj >  250.0*GeV)  _h["zjj-ew"]->fillBin(0);
        if (zpt > 20.0*GeV && !ngapjets && ptbalance2 < 0.15 && mjj > 1000.0*GeV)  _h["zjj-ew"]->fillBin(1);
      }
      else {
        _h["zjj"]->fillBin(0);
        if (jet1pt > 85.0*GeV && jet2pt > 75.0*GeV)  _h["zjj"]->fillBin(1);
        if (zpt > 20.0*GeV && ngapjets == 0 && ptbalance2 < 0.15 && mjj > 250.0*GeV)  _h["zjj"]->fillBin(2);
        if (zpt > 20.0*GeV && ngapjets && ptbalance3 < 0.15 && mjj > 250.0*GeV)  _h["zjj"]->fillBin(3);
        if (mjj > 1000.0*GeV)  _h["zjj"]->fillBin(4);
        if (zpt > 20.0*GeV && !ngapjets && ptbalance2 < 0.15 && mjj > 1000.0*GeV)  _h["zjj"]->fillBin(2);
      }
     
    }
    
    
    
     /// Normalise histograms etc., after the run
    void finalize() {

      double factor = crossSection()/femtobarn/sumOfWeights();
      scale(_h, factor);
    }
    
    bool _isBetween(const Jet probe, const Jet boundary1, const Jet boundary2) {
      double y_p = probe.rapidity();
      double y_b1 = boundary1.rapidity();
      double y_b2 = boundary2.rapidity();

      double y_min = std::min(y_b1, y_b2);
      double y_max = std::max(y_b1, y_b2);

      if (y_p > y_min && y_p < y_max) return true;
      else return false;
    }

    ///@}

  private:

    size_t _mode;

    /// @name Histograms
    ///@{
     map<string, Histo1DPtr> _h;
    ///@}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2017_I1627873);

}

