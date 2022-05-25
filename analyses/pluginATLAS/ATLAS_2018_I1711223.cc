// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"


namespace Rivet {


  /// @brief Electroweak WZjj production cross section at 13 TeV
  class ATLAS_2018_I1711223 : public Analysis {
  public:
   
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1711223);

    /// @name Analysis methods
    //@{
    
    /// Book histograms and initialise projections before the run
    void init() {

      // Get photons to dress leptons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Electrons and muons in Fiducial PS
      PromptFinalState leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      leptons.acceptTauDecays(false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, Cuts::open(), true);                                       // useDecayPhotons=true -- useJetClustering? auto-set to false?
      declare(dressedleptons, "DressedLeptons");

      // Prompt neutrinos (yikes!)
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(false);
      declare(neutrinos, "Neutrinos");
      MSG_WARNING("\033[91;1mLIMITED VALIDITY - check info file for details!\033[m");

      //Jets
    
    	// Muons
    	PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true); // true = use muons from prompt tau decays
    	DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, Cuts::abseta < 5.0, true);

    	// Electrons
    	PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true); // true = use electrons from prompt tau decays
    	DressedLeptons all_dressed_el(photons, bare_el, 0.1, Cuts::abseta < 5.0, true);

    	//Jet forming
    	VetoedFinalState vfs(FinalState(Cuts::abseta < 5));
    	vfs.addVetoOnThisFinalState(all_dressed_el);
    	vfs.addVetoOnThisFinalState(all_dressed_mu);
          
    	FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
    	declare(jets, "Jets");

      // Book auxiliary histograms
      book(_h["MTWZ"],         "_mTWZ", refData( 6, 1, 1));
      book(_h["sumpt"],       "_sumpT", refData( 8, 1, 1));
      book(_h["dphiWZ"],     "_dphiWZ", refData(10, 1, 1));
      book(_h["Njets_VBS"],   "_njets", refData(12, 1, 1));
      book(_h["mjj"],           "_mjj", refData(14, 1, 1));
      book(_h["dyjj"],       "_dRapjj", refData(16, 1, 1));
      book(_h["dphijj"],     "_dPhijj", refData(18, 1, 1));
      book(_h["Njets_gap"], "_gapJets", refData(20, 1, 1));

      // book output bar charts
      book(_s["MTWZ"],       6, 1, 1);
      book(_s["sumpt"],      8, 1, 1);
      book(_s["dphiWZ"],    10, 1, 1);
      book(_s["Njets_VBS"], 12, 1, 1);
      book(_s["mjj"],       14, 1, 1);
      book(_s["dyjj"],      16, 1, 1);
      book(_s["dphijj"],    18, 1, 1);
      book(_s["Njets_gap"], 20, 1, 1);
      
    }


    void analyze(const Event& event) {

      const Particles& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.5);

      int i, j, k;
      double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
      double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
      double WeightZ1, WeightZ2, WeightZ3;
      double WeightW1, WeightW2, WeightW3;
      double M1, M2, M3;
      double WeightTotal1, WeightTotal2, WeightTotal3;

      //---Fiducial PS: assign leptons to W and Z bosons using Resonant shape algorithm
      if (dressedleptons.size() < 3 || neutrinos.size() < 1) vetoEvent;                                                              
 
      //--- count num of electrons and muons
      int Nel = 0, Nmu = 0;
      for (const Particle& l : dressedleptons) {
        if (l.abspid() == 11)  ++Nel;
        if (l.abspid() == 13)  ++Nmu;
      }

      int icomb=0;
      // try Z pair of leptons 01                                                              
    if ( (dressedleptons[0].pid() ==-(dressedleptons[1].pid()))  && (dressedleptons[2].pid()*neutrinos[0].pid()< 0) && (dressedleptons[2].abspid()==neutrinos[0].abspid()-1)) {
        MassZ01 = (dressedleptons[0].momentum() + dressedleptons[1].momentum()).mass();
        MassW2 = (dressedleptons[2].momentum() + neutrinos[0].momentum()).mass();
        icomb = 1;
      }

      // try Z pair of leptons 02
      if ( (dressedleptons[0].pid()==-(dressedleptons[2].pid()))  && (dressedleptons[1].pid()*neutrinos[0].pid()< 0) && (dressedleptons[1].abspid()==neutrinos[0].abspid()-1)) {
        MassZ02 = (dressedleptons[0].momentum() + dressedleptons[2].momentum()).mass();
        MassW1 = (dressedleptons[1].momentum() + neutrinos[0].momentum()).mass();
        icomb = 2;
      }
      // try Z pair of leptons 12
      if ( (dressedleptons[1].pid()==-(dressedleptons[2].pid())) && (dressedleptons[0].pid()*neutrinos[0].pid()< 0) && (dressedleptons[0].abspid()==neutrinos[0].abspid()-1)) {
        MassZ12 = (dressedleptons[1].momentum() + dressedleptons[2].momentum()).mass();
        MassW0 = (dressedleptons[0].momentum() + neutrinos[0].momentum()).mass();
        icomb = 3;
      }
 
      if (icomb<=0)  vetoEvent;


      WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ) {
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ) {
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ) {
        i = 1; j = 2; k = 0;
      }

      FourMomentum Zlepton1 = dressedleptons[i].momentum();
      FourMomentum Zlepton2 = dressedleptons[j].momentum();
      FourMomentum Wlepton  = dressedleptons[k].momentum();
      FourMomentum Zboson   = dressedleptons[i].momentum()+dressedleptons[j].momentum();
      FourMomentum Wboson   = dressedleptons[k].momentum()+neutrinos[0].momentum();

      double cosLepNeut;
      double Wboson_mT = 0.;
      double norm = Wlepton.pT() * neutrinos[0].pt();
      if(norm != 0 ) {
        cosLepNeut = ( Wlepton.px()*neutrinos[0].px() + Wlepton.py()*neutrinos[0].py() )/norm ;
        if (1-cosLepNeut >= 0. ) Wboson_mT = sqrt( 2 * Wlepton.pT() * neutrinos[0].pt() * (1-cosLepNeut ) );
      }

      //---- CUTS (based on Table 1 WZ: 36.1 fb-1)----//
      if (Wlepton.pT() <= 20*GeV || Zlepton1.pT() <= 15*GeV || Zlepton2.pT() <= 15*GeV)     vetoEvent;      
      if (Wlepton.abseta() >= 2.5 || Zlepton1.abseta() >= 2.5 || Zlepton2.abseta() >= 2.5)  vetoEvent;
      if (fabs(Zboson.mass()/GeV - MZ_PDG) >= 10.) vetoEvent;
      if (Wboson_mT <= 30*GeV)                     vetoEvent;
      if (deltaR(Zlepton1, Zlepton2) <= 0.2)       vetoEvent;
      if (deltaR(Zlepton1, Wlepton)  <= 0.3)       vetoEvent;
      if (deltaR(Zlepton2, Wlepton)  <= 0.3)       vetoEvent;

      double WZ_pt = (Zlepton1.pt() + Zlepton2.pt() + Wlepton.pt() + neutrinos[0].pt())/GeV;
      double WZ_px = (Zlepton1.px() + Zlepton2.px() + Wlepton.px() + neutrinos[0].px())/GeV;
      double WZ_py = (Zlepton1.py() + Zlepton2.py() + Wlepton.py() + neutrinos[0].py())/GeV;
      double mTWZ = sqrt( pow(WZ_pt, 2) - ( pow(WZ_px, 2) + pow(WZ_py,2) ) );
      double sumptleptons = (Zlepton1.pt() + Zlepton2.pt() + Wlepton.pt())/GeV;
      double dPhiWZTruth = acos(cos(Zboson.phi()-Wboson.phi()));

    
      
      //---- Jet CUTS----//
      ifilter_discard(jets, [&](const Jet& j) {
        return deltaR(j, Zlepton1) < 0.3 || deltaR(j, Zlepton2) < 0.3 || deltaR(j, Wlepton) < 0.3;
      });
      if (jets.size() < 2)  vetoEvent;  
      if (jets[0].pT() < 40*GeV)  vetoEvent;  

      // Selection of the second jet as the second highest pT jet and in opposite hemisphere with the fisrt jet
      FourMomentum jet_lead = jets[0].mom();
      FourMomentum jet_sublead;
      bool foundVBSJetPair = false;
      for (const Jet& jet : jets) {
        if(jet.pT() > 40*GeV && jet.eta()*jets[0].eta() < 0.) {
          jet_sublead = jet.mom();
          foundVBSJetPair = true;
          break;
        }
      }
      if (!foundVBSJetPair)  vetoEvent;
     
      const double mJJ = (jet_lead + jet_sublead).mass()/GeV;
      const double dphi_jj = acos(cos(jet_lead.phi() - jet_sublead.phi()));
      const double dyjj = fabs(jet_lead.rap() - jet_sublead.rap());

      //Plots in the SR
      if (mJJ < 500*GeV) vetoEvent;

      const size_t njets40 = filter_select(jets, Cuts::pT > 40*GeV).size();
      fillWithOverflow("Njets_VBS", njets40, 5.1);

      const double y_min = std::min(jet_lead.rap(), jet_sublead.rap());
      const double y_max = std::max(jet_lead.rap(), jet_sublead.rap());
      const size_t njetsGap = count(jets, [&](const Jet& j) { 
        return  (j.rap() > y_min && j.rap() < y_max);
      });
      fillWithOverflow("Njets_gap", njetsGap, 3.1);

      fillWithOverflow("MTWZ", mTWZ, 551);
      fillWithOverflow("sumpt", sumptleptons, 501);
      fillWithOverflow("mjj", mJJ, 2001);

      _h["dphiWZ"]->fill(dPhiWZTruth);
      _h["dyjj"]->fill(dyjj);
      _h["dphijj"]->fill(dphi_jj);

    }


    void fillWithOverflow(const string& tag, const double value, const double overflow) {
      if (value < overflow)  _h[tag]->fill(value);
      else                   _h[tag]->fill(overflow);
    }


    void finalize() {

      scale(_h, crossSectionPerEvent() / femtobarn);
      // unfortunately, no differential cross-sections were measured in this analysis
      for (auto &item : _h)  barchart(item.second, _s[item.first]);

    }


 //@}

  private:


    /// @name Histograms
    //@{

    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;

    //@}

    double MZ_PDG = 91.1876;
    double MW_PDG = 80.385;
    double GammaZ_PDG = 2.4952;
    double GammaW_PDG = 2.085;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1711223);
}
