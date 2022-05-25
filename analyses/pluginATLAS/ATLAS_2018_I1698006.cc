#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"



namespace Rivet {
  
  /// @brief ATLAS pTmiss+gamma measurement at 13 TeV 
  class ATLAS_2018_I1698006 : public Analysis {
  public:
    
    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1698006);
    
    /// @name Analysis methods
    //@{
    void init() {

      // Get options
      // Default (OFF) uses the extended phase space of the paper.
      // LVETO ON will add an additonal lepton veto to reject non-Zgamma events
      // (for conservative BSM limits, for example).
      _mode = 0;
      if ( getOption("LVETO") == "ON" ) _mode = 1;

      //prompt photons
      const Cut photoncut = Cuts::abspid == PID::PHOTON && Cuts::Et > 150*GeV && Cuts::abseta < 2.37;
      const PromptFinalState photon_fs(photoncut);
      declare(photon_fs, "Photons");

      //missing energy (prompt neutrinos)
      declare(InvisibleFinalState(true), "MET");
      
      if (_mode==1) {
	FinalState allLeps(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
	FinalState photons(Cuts::abspid == PID::PHOTON);
	PromptFinalState promptLeps(allLeps);
	Cut dressedLep_cuts = (Cuts::abseta < 2.7) && (Cuts::pT > 7*GeV);
	const DressedLeptons dressedLeps(photons, promptLeps, 0.1, dressedLep_cuts, true);
	declare(dressedLeps, "dressedLeptons");
      }
      
      //jets. run the jet finder on a final state without the prompt photons, and without neutrinos or muons
      VetoedFinalState jet_fs(Cuts::abseta > 4.5);
      jet_fs.addVetoOnThisFinalState(photon_fs);
      FastJets fastjets(jet_fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(fastjets, "Jets");


      //books histograms
      //fig.4a
      book(_h["Et_inc"],2,1,1);
      //fig.4b
      book(_h["Et_exc"],3,1,1);
      //fig.5a
      book(_h["pT_inc"],4,1,1);
      //fig.5b
      book(_h["pT_exc"],5,1,1);
      //fig.6
      book(_h["Njets"],6,1,1);

    }
    
    void analyze(const Event& event) {


      const Particles& photons = apply<PromptFinalState>(event,"Photons").particlesByPt();
      const Jets& jets = apply<FastJets>(event,"Jets").jetsByPt(Cuts::pT > 50*GeV);
      const FinalState& metfs = apply<InvisibleFinalState>(event,"MET");
      Vector3 met_vec;
      for (const Particle& p : metfs.particles()) met_vec += p.mom().perpVec();


      if (_mode==1) {
	const vector<DressedLepton> &dressedLeptons = apply<DressedLeptons>(event, "dressedLeptons").dressedLeptons();
	if (dressedLeptons.size() > 0) vetoEvent;
      }
      
      
      //Nγ==1 and Emiss > 150 GeV
      if (met_vec.mod() > 150*GeV && photons.size()==1){

	  //inclusive case (Njet>=0)
	  if (jets.size()>=0){
	    bool dR_veto = any(jets, DeltaRLess(photons[0], 0.3));
	    if (not dR_veto) {
	      double Et_photon = photons[0].Et()/GeV;
	      _h["Et_inc"]->fill(Et_photon);
	      //fill in missing energy (neutrino) pT histogram (inclusive)
	      _h["pT_inc"]->fill(met_vec.mod()/GeV);
	    }
	  }

	  //exclusive case (Njet==0)
	  if (jets.size() == 0){
	    double Et_photon = photons[0].Et()/GeV;
	    _h["Et_exc"]->fill(Et_photon);
	     //fill in missing energy (neutrino) pT histogram (exclusive)
	    _h["pT_exc"]->fill(met_vec.mod()/GeV);
	  }

	  _h["Njets"]->fill(jets.size());
	  
      }
      
    }

    void finalize() {
      
      const double sf = crossSection()/femtobarn/sumOfWeights();
      scale(_h, sf);

    }

    //@}

    private:
    map<string, Histo1DPtr> _h;
    size_t _mode;
    
  };

  // Magic required by the plugin system 
  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1698006);
  
}
