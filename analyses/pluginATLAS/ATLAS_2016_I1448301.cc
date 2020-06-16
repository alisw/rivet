// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {


  /// @brief Z/gamma cross section measurement at 8 TeV
  class ATLAS_2016_I1448301 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1448301);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // Default tries to fill everything
      // NU fills only the MET plots
      // EL fills only the Electron plots and combined plots assuming lepton univerality
      // MU fills only the Muon plots and combined plots assuming lepton univerality
      // LL fills electron and Muon plots and combined plots from correct average
      _mode = 0;
      if ( getOption("LMODE") == "NU" ) _mode = 1;
      if ( getOption("LMODE") == "EL" ) _mode = 2;
      if ( getOption("LMODE") == "MU" ) _mode = 3;
      if ( getOption("LMODE") == "LL" ) _mode = 4;

      // Prompt photons
      const Cut photoncut = Cuts::abspid == PID::PHOTON && Cuts::pT > 15*GeV && Cuts::abseta < 2.37;
      PromptFinalState photon_fs(photoncut);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const IdentifiedFinalState allphoton_fs(PID::PHOTON); // photons used for lepton dressing
      const Cut leptoncut = Cuts::pT > 25*GeV && Cuts::abseta < 2.47;
      const DressedLeptons dressedelectron_fs(allphoton_fs, bareelectron_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing
      const DressedLeptons dressedmuon_fs(allphoton_fs, baremuon_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");

      // MET (prompt neutrinos)
      VetoedFinalState ivfs;
      ivfs.addVetoOnThisFinalState(VisibleFinalState());
      declare(PromptFinalState(ivfs), "MET");

      // Jets
      VetoedFinalState jet_fs;
      jet_fs.vetoNeutrinos();
      jet_fs.addVetoPairId(PID::MUON);
      const FastJets fastjets(jet_fs, FastJets::ANTIKT, 0.4);
      declare(fastjets, "Jets");

      // Histograms

      // MET
      if (_mode == 0 || _mode == 1){
        book(_h["vvg"],     2, 1, 1);
        book(_h["vvgg"],    4, 1, 1);
        book(_h["pT"],      7, 1, 1);
        book(_h["pT_0jet"], 8, 1, 1);
      }
	
      // always book e and mu in charged lepton modes; there are sometimes 4 leptons.
      if (_mode != 1){
	// electron
        book(_h["eeg"],  1, 1, 1);
        book(_h["eegg"], 3, 1, 1);
        // muon
	book(_h["mmg"],  1, 1, 2);
	book(_h["mmgg"], 3, 1, 2);

        // combined
        book(_h["llgg"], 3, 1, 3);
        book(_h["llg"], 1, 1, 3);
        book(_h["pT"], 5, 1, 1);
        book(_h["pT_0jet"], 6, 1, 1);
        book(_h["M"], 9, 1, 1);
        book(_h["M_0jet"], 10, 1, 1);
        book(_h["Njets"], 11, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get objects
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "Electrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "Muons").dressedLeptons();

      const Particles& photons = apply<PromptFinalState>(event, "Photons").particlesByPt();
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt();

      const FinalState& metfs = apply<PromptFinalState>(event, "MET");
      Vector3 met_vec;
      for (const Particle& p : metfs.particles()) met_vec += p.mom().perpVec();
      
      if (met_vec.mod() >= 100*GeV && !photons.empty() && _mode < 2){

	if (photons.size() > 1) { // nu nu y y
	
	  bool yy_veto = false;
	  yy_veto |= photons[0].pT() < 22*GeV;
	  yy_veto |= photons[1].pT() < 22*GeV;
	  yy_veto |= met_vec.mod() < 110*GeV;
	  const double yyPhi = (photons[0].momentum() + photons[1].momentum()).phi();
	  yy_veto |= fabs(yyPhi - met_vec.phi()) < 2.62 || fabs(yyPhi - met_vec.phi()) > 3.66;
	  yy_veto |= deltaR(photons[0], photons[1]) < 0.4;
	  
	  // Photon isolation calculated by jets, count jets
	  Jet ph0_jet, ph1_jet;
	  double min_dR_ph0_jet = 999., min_dR_ph1_jet = 999.;
	  size_t njets = 0;
	  for (const Jet& j : jets) {
	    if (j.pT() > 30*GeV && j.abseta() < 4.5) {
	      if (deltaR(j, photons[0]) > 0.3 && deltaR(j, photons[1]) > 0.3)  ++njets;
	    }
	    if (deltaR(j, photons[0]) < min_dR_ph0_jet) {
	      min_dR_ph0_jet = deltaR(j, photons[0]);
	      ph0_jet = j;
	    }
	    if (deltaR(j, photons[1]) < min_dR_ph1_jet) {
	      min_dR_ph1_jet = deltaR(j, photons[1]);
	      ph1_jet = j;
	    }
	  }
	  double photon0iso = 0., photon1iso = 0.;
	  if (min_dR_ph0_jet < 0.4)  photon0iso = ph0_jet.pT() - photons[0].pT();
	  if (min_dR_ph1_jet < 0.4)  photon1iso = ph1_jet.pT() - photons[1].pT();
	  yy_veto |= photon0iso/photons[0].pT() > 0.5;
	  yy_veto |= photon1iso/photons[1].pT() > 0.5;
	  
	  if (!yy_veto) {
	    _h["vvgg"]->fill(0.5);
	    if (!njets)  _h["vvgg"]->fill(1.5);
	  }
	} // end of nu nu y y section
      
      
	if ((photons[0].pT() >= 130*GeV)  &&
	    (fabs(fabs(deltaPhi(photons[0], met_vec)) - 3.14) <= 1.57)) {
      
	  // Photon isolation calculated by jets, count jets
	  Jet ph_jet;
	  double min_dR_ph_jet = 999.;
	  size_t njets = 0;
	  for (const Jet& j : jets) {
	    if (j.pT() > 30*GeV && j.abseta() < 4.5) {
	      if (deltaR(j, photons[0]) > 0.3)  ++njets;
	    }
	    if (deltaR(j, photons[0]) < min_dR_ph_jet) {
	      min_dR_ph_jet = deltaR(j, photons[0]);
	      ph_jet = j;
	    }
	  }
	  double photoniso = 0;
	  if (min_dR_ph_jet < 0.4)  photoniso = ph_jet.pT() - photons[0].pT();
	  if (photoniso/photons[0].pT() > 0.5)  vetoEvent;
	  
	  const double pTgamma = photons[0].pT()/GeV;
	  _h["pT"]->fill(pTgamma);
	  _h["vvg"]->fill(0.5);
	  if (!njets) {
	    _h["vvg"]->fill(1.5);
	    _h["pT_0jet"]->fill(pTgamma);
	  }
	  	  
	}
      }  // end of nu nu y (y) section

      // Dilepton candidate
      bool el = false;
      if ( (_mode != 1) &&
	   (( electrons.size() >= 2 && _mode != 3 ) || 
	    ( muons.size()     >= 2 && _mode != 2 ) )) {

	vector<DressedLepton> lep_p, lep_m;

	// Sort the dressed leptons by pt
	if (electrons.size() >= 2) {
	  el = true;
	  sortByPt(electrons);
	  for (const DressedLepton& lep : electrons) {
	    if (lep.charge() > 0.)  lep_p.push_back(lep);
	    if (lep.charge() < 0.)  lep_m.push_back(lep);
	  }	
	} else {
	  sortByPt(muons);
	  for (const DressedLepton& lep : muons) {
	    if (lep.charge() > 0.)  lep_p.push_back(lep);
	    if (lep.charge() < 0.)  lep_m.push_back(lep);
	  }	
	}
      

	if (!lep_p.empty() && !lep_m.empty() &&
	    (lep_p[0].abspid() == lep_m[0].abspid()) &&
	    ((lep_p[0].momentum() + lep_m[0].momentum()).mass() >= 40*GeV)){

	  // Photon lepton overlap removal
	  if (photons.empty())  vetoEvent;
	  
	  if (photons.size() > 1) {
	    
	    bool veto = false;
	    veto |= deltaR(photons[0], lep_p[0]) < 0.4;
	    veto |= deltaR(photons[0], lep_m[0]) < 0.4;
	    veto |= deltaR(photons[1], lep_p[0]) < 0.4;
	    veto |= deltaR(photons[1], lep_m[0]) < 0.4;
	    veto |= deltaR(photons[0], photons[1]) < 0.4;
	    
	    Jet ph0_jet, ph1_jet;
	    double min_dR_ph0_jet = 999., min_dR_ph1_jet=999.;
	    int njets = 0;
	    for (const Jet& j : jets){
	      if (j.pT() > 30*GeV && j.abseta() < 4.5) {
		if (deltaR(j, lep_p[0]) > 0.3 && deltaR(j, lep_m[0]) > 0.3) {
		  if (deltaR(j, photons[0]) > 0.3 && deltaR(j, photons[1]) > 0.3 )  ++njets;
		}
	      }
	      if (deltaR(j, photons[0]) < min_dR_ph0_jet) {
		min_dR_ph0_jet = deltaR(j, photons[0]);
		ph0_jet = j;
	      }
	      if (deltaR(j, photons[1]) < min_dR_ph1_jet) {
		min_dR_ph1_jet = deltaR(j, photons[1]);
		ph1_jet = j;
	      }
	    }
	    double photon0iso = 0, photon1iso = 0;
	    if (min_dR_ph0_jet < 0.4) photon0iso = ph0_jet.pT() - photons[0].pT();
	    if (min_dR_ph1_jet < 0.4) photon1iso = ph1_jet.pT() - photons[1].pT();
	    veto |= photon0iso/photons[0].pT() > 0.5;
	    veto |= photon1iso/photons[1].pT() > 0.5;
	    
	    // Fill plots
	    // ee and mm need doing.
	    if (!veto) {
	      _h["llgg"]->fill(0.5);
	      if (el) {
		_h["eegg"]->fill(0.5);
	      } else {
		_h["mmgg"]->fill(0.5);
	      }

	      if (!njets) {
		_h["llgg"]->fill(1.5);
		if (el) {
		  _h["eegg"]->fill(1.5);
		} else {
		  _h["mmgg"]->fill(1.5);
		}
	      }
	    }
	  }
	  
	  if (deltaR(photons[0], lep_p[0]) < 0.7)  vetoEvent;
	  if (deltaR(photons[0], lep_m[0]) < 0.7)  vetoEvent;
	  
	  // Photon isolation calculated by jets, count jets
	  Jet ph_jet;
	  double min_dR_ph_jet = 999.;
	  size_t njets = 0;
	  for (const Jet& j : jets) {
	    if (j.pT() > 30*GeV && j.abseta() < 4.5) {
	      if (deltaR(j, lep_p[0]) > 0.3 && deltaR(j, lep_m[0]) > 0.3 && deltaR(j, photons[0]) > 0.3)  ++njets;
	    }
	    if (deltaR(j, photons[0]) < min_dR_ph_jet) {
	      min_dR_ph_jet = deltaR(j, photons[0]);
	      ph_jet = j;
	    }
	  }
	  
	  double photoniso = 0;
	  if (min_dR_ph_jet < 0.4)  photoniso = ph_jet.pT() - photons[0].pT();
	  if (photoniso/photons[0].pT() > 0.5)  vetoEvent;
	  
	  
	  // Fill plots
	  const double pTgamma = photons[0].pT()/GeV;
	  const double mllgamma = (lep_p[0].momentum() + lep_m[0].momentum() + photons[0].momentum()).mass()/GeV;
	  
	  _h["pT"]->fill(pTgamma);
	  _h["M"]->fill(mllgamma);
	  _h["Njets"]->fill(njets < 3? njets : 3);
	  
	  _h["llg"]->fill(0.5);
	  if (el) {
	    _h["eeg"]->fill(0.5);
	  } else {
	    _h["eeg"]->fill(0.5);
	  }

	  if (!njets) {
	    _h["pT_0jet"]->fill(pTgamma);
	    _h["M_0jet"]->fill(mllgamma);
	    _h["llg"]->fill(1.5);
	    if (el) {
	      _h["eeg"]->fill(1.5);
	    } else {
	      _h["mmg"]->fill(1.5);
	    }
	  }
	}
      }
    } // end of analysis


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/femtobarn/sumOfWeights();
      for (const auto& kv : _h) scale(kv.second, sf);
      // if we are running both e and mu, the combined lepton histos
      // need to be divided by two to get the average
      if (_mode == 0 || _mode == 4){
        scale(_h["llgg"], 0.5);
        scale(_h["llg"], 0.5);
        scale(_h["pT"], 0.5);
        scale(_h["pT_0jet"], 0.5);
        scale(_h["M"], 0.5);
        scale(_h["M_0jet"], 0.5);
        scale(_h["Njets"], 0.5);
      }
    }

    //@}


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode;

  private:

    /// Histograms
    map<string, Histo1DPtr> _h;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1448301);

}
