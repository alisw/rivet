// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief WWW cross-section at 8 TeV, 3L mode
  class ATLAS_2016_I1492320 : public Analysis {
  public:

    // Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1492320);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "3L" ) _mode = 1;
      if ( getOption("LMODE") == "2L2J" ) _mode = 2;


      // Charged leptons within acceptance
      const PromptFinalState chLep_fid = PromptFinalState(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      const PromptFinalState photon_fs = PromptFinalState(Cuts::abspid == PID::PHOTON);
      const DressedLeptons dressed_leps(photon_fs, chLep_fid, 0.1, Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      declare(dressed_leps, "DressedLeptons");

      const DressedLeptons dressed_leps2(photon_fs, chLep_fid, 0.1, Cuts::pT > 10*GeV);
      declare(dressed_leps2, "DressedLeptons2");


      // Jets, anti-kt 0.4
      VetoedFinalState fsJets(FinalState(Cuts::abseta < 7.0)); //final state for jet finding: veto leptons and neutrinos
      fsJets.vetoNeutrinos();
      fsJets.addVetoOnThisFinalState(photon_fs);
      fsJets.addVetoOnThisFinalState(chLep_fid);
      declare(FastJets(fsJets, FastJets::ANTIKT, 0.4), "Jets");

      // b hadrons for b-tagging
      declare(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "Bhadrons");

      // Missing momentum
      declare(MissingMomentum(), "MET");

      // Histograms
      if (_mode != 2){
	_h_fiducial_3l = bookCounter("d01-x01-y01");
      }
      if (_mode != 1){
	_h_2l2j = bookCounter("d01-x01-y02");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the dressed leptons, sorted by pT of their constituent bare lepton (!!)
      vector<DressedLepton> _vbs_lep = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      if (_vbs_lep.size() == 3 && _mode != 2) {
	std::sort(_vbs_lep.begin(), _vbs_lep.end(), [](const DressedLepton& l1, const DressedLepton& l2) {
	    return (l1.constituentLepton().pT() > l2.constituentLepton().pT());
	  });
	
	// Get the jets
	const Jets& _vbs_jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.5);
	if (_vbs_jets.size() <= 1) {
	
	  // Determine nsfos pairs for channel classification
	  int nSFOS = 0;
	  for (size_t i = 0; i < _vbs_lep.size(); ++i) {
	    const double ch_l0 = _vbs_lep[i].charge();
	    for (size_t j = i + 1; j < _vbs_lep.size(); ++j) {
	      const double ch_l1 = _vbs_lep[j].charge();
	      if (_vbs_lep[i].abspid() == _vbs_lep[j].abspid() && ch_l0*ch_l1 < 0) ++nSFOS;
	    }
	  }
	  
	  double minDRll = DBL_MAX, mSFOS_MinDiff = DBL_MAX, meeSS_MinDiff = DBL_MAX, mSF_min = DBL_MAX;
	  for (size_t i = 0; i < _vbs_lep.size(); ++i) {
	    const double ch_l0 = _vbs_lep[i].charge();
	    for (size_t j = i + 1; j < _vbs_lep.size(); ++j) {
	      const double ch_l1 = _vbs_lep[j].charge();
	      const bool samesign = ch_l0*ch_l1 > 0;
	      
	      // Update min dR between leptons
	      minDRll = min(minDRll, deltaR(_vbs_lep[i], _vbs_lep[j]));
	      
	      // Require same flavour
	      if (_vbs_lep[i].abspid() != _vbs_lep[j].abspid()) continue;
	      
	      // SF dilepton mass (used several times)
	      const double mSF = (_vbs_lep[i].momentum() + _vbs_lep[j].momentum()).mass();
	      
	      // Determine min for all same-flavor pairs
	      mSF_min = min(mSF, mSF_min);
	      
	      // Determine min for all m_ee same-sign pairs
	      if (_vbs_lep[i].abspid() == PID::ELECTRON && samesign) {
		if (fabs(mSF-ZMASS) < fabs(meeSS_MinDiff-ZMASS)) meeSS_MinDiff = mSF;
	      }
	      
	      // Determine min for all mSFOS pairs
	      if (!samesign && fabs(mSF-ZMASS) < abs(mSFOS_MinDiff-ZMASS)) mSFOS_MinDiff = mSF;
	    }
	  }
	  
	  bool setVeto = false;
	  if (minDRll < 0.1) setVeto = true;
	  if (nSFOS == 0 && mSF_min < 20*GeV) setVeto = true;
	  if (nSFOS == 0 && fabs(meeSS_MinDiff - ZMASS) < 15*GeV) setVeto = true;
	  if (nSFOS == 1 && ((ZMASS - mSFOS_MinDiff) < 35*GeV && (mSFOS_MinDiff - ZMASS) < 20*GeV)) setVeto = true;
	  if (nSFOS == 2 && fabs(mSFOS_MinDiff - ZMASS) < 20*GeV) setVeto = true;
	  
	  if (!setVeto) {
	    const Vector3& met = -1.0 * apply<MissingMomentum>(event, "MET").vectorEt();
	    if (nSFOS == 1 && met.mod() < 45*GeV) setVeto = true;
	    if (nSFOS == 2 && met.mod() < 55*GeV) setVeto = true;
	    
	    if (!setVeto) {
	      const double dPhi = deltaPhi((_vbs_lep[0].momentum() + _vbs_lep[1].momentum() + _vbs_lep[2].momentum()), met);
	      if (dPhi < 2.5) setVeto = true;
	    }
	    // Fill histo
	    if (!setVeto) {
	      _h_fiducial_3l->fill(event.weight());
	    }
	  }
	}
      }

      if (_mode != 1){
	// Get leptons
	vector<DressedLepton> leps = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
	if (leps.size() >= 2) {
	  // Sort the dressed leptons by pt of their constituent lepton (bare lepton pt)
	  std::sort(leps.begin(), leps.end() ,
		    [](const DressedLepton& l1, const DressedLepton& l2) {
		      return (l1.constituentLepton().pT() > l2.constituentLepton().pT()); });
	  if (leps[0].pT() < 30*GeV || leps[0].abseta() > 2.5)  vetoEvent;
	  if (leps[1].pT() < 30*GeV || leps[1].abseta() > 2.5)  vetoEvent;
	  
	  // Get jets
	  const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 15*GeV);
	  
	  // Find min dilepton DR and mass
	  double minDRll = DBL_MAX, mll = DBL_MAX;
	  for (size_t i = 0; i < leps.size(); ++i) {
	    for (size_t j = i + 1; j < leps.size(); ++j) {
	      minDRll = min(minDRll, deltaR(leps[i], leps[j]));
	      mll = min(mll, (leps[i].momentum() + leps[j].momentum()).mass());
	    }
	  }
	  if (minDRll < 0.1) vetoEvent;
	  if (mll < 40*GeV) vetoEvent;
	  
	  // Require same-sign leading leptons
	  if (leps[0].charge()*leps[0].charge() < 0) vetoEvent;
	  
	  // Veto di-electron combinations within 10 GeV of the Z mass
	  if (fabs(mll - 91.188*GeV) < 10*GeV && leps[0].abspid() == PID::ELECTRON && leps[1].abspid() == PID::ELECTRON) vetoEvent;
	  
	  // Now jet cuts
	  if (jets.size() < 2) vetoEvent;
	  if (jets[0].pT() < 30*GeV || jets[0].abseta() > 2.5) vetoEvent;
	  if (jets[1].pT() < 20*GeV || jets[1].abseta() > 2.5) vetoEvent;
	  
	  // Find closest jet/lepton pair and veto if too close in phi or too far in eta
	  double minDRLepJets = DBL_MAX;
	  for (const Jet& jet : jets) {
	    for (const Particle& lep : leps) minDRLepJets = min(minDRLepJets, deltaR(lep, jet));
	  }
	  if (minDRLepJets < 0.3) vetoEvent;
	  if (fabs(deltaEta(jets[0], jets[1])) > 1.5) vetoEvent;
	  
	  // Dijet mass requirement
	  double mjj = (jets[0].momentum() + jets[1].momentum()).mass();
	  if (mjj < 65 || mjj > 105)  vetoEvent;
	  if (!inRange(mjj, 65*GeV, 105*GeV)) vetoEvent;
	  
	  // Veto if any good jets are b-jets
	  const Particles& bhadrons = apply<HeavyHadrons>(event, "Bhadrons").bHadrons();
	  for (const Jet& j : jets) {
	    if (j.abseta() > 2.5) continue; // outside acceptance of b-tagging
	    const bool isbJet = any(bhadrons, deltaRLess(j, 0.3));
	    if (isbJet) vetoEvent;
	  }
	  
	  // MET vetoing for non-muon events
	  const MissingMomentum& met = apply<MissingMomentum>(event, "MET");
	  if (met.vectorEt().mod() < 55*GeV && (leps[0].abspid() != PID::MUON || leps[1].abspid() != PID::MUON)) vetoEvent;

	  // Fill counter
	  _h_2l2j->fill(event.weight());
	}
      }
    }



    /// Normalise histograms etc., after the run
    void finalize() {
      if (_mode != 2){
	scale(_h_fiducial_3l, crossSection()/sumOfWeights()/femtobarn);
      }
      if (_mode != 1){
	scale(_h_2l2j, crossSection()/sumOfWeights()/femtobarn);
      }
    }

    //@}

  protected:

    size_t _mode;


  private:

    /// @name Histograms
    //@{
    const double ZMASS = 91.1876*GeV;
    CounterPtr _h_fiducial_3l;
    CounterPtr _h_2l2j;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1492320);

}
