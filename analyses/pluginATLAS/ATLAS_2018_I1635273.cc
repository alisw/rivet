#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

  /// @brief W + jets production at 8 TeV
  class ATLAS_2018_I1635273 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1635273);

	  // Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // Default uses electrons
      // EL looks for W->ev candidate
      // MU looks for W->mv candidate
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;

      FinalState fs;
      Cut cuts = Cuts::pT > 25*GeV && Cuts::abseta < 2.5;

      // Get photons to dress leptons
      // (Paper says "radiated photons", but there was
      // no promptness requirement in the analysis code)
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Get dressed leptons
      PromptFinalState leptons(Cuts::abspid == (_mode? PID::MUON : PID::ELECTRON), false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, cuts, true);
      declare(dressedleptons, "DressedLeptons");

      // Get neutrinos for MET calculation
      declare(InvisibleFinalState(true), "InvFS"); // true = only allow prompt invisibles

      // jets
      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE);
      declare(jets, "Jets");

      // book histograms
      book(_h["N_incl_pb"],            1, 1, 1);
      book(_h["HT_1j_fb"],             6, 1, 1);
      book(_h["W_pt_1j_fb"],          11, 1, 1);
      book(_h["jet_pt1_1j_fb"],       16, 1, 1); 
      book(_h["jet_y1_1j_fb"],        21, 1, 1);
      book(_h["jet_pt2_2j_fb"],       26, 1, 1);
      book(_h["jet_y2_2j_fb"],        28, 1, 1);
      book(_h["DeltaRj12_2j_fb"],     30, 1, 1);
      book(_h["jet_mass12_2j_fb"],    32, 1, 1);
      book(_h["N_pb"],                34, 1, 1);
      book(_h["HT_2j_fb"],            36, 1, 1);
      book(_h["W_pt_2j_fb"],          41, 1, 1);
      book(_h["jet_pt1_2j_fb"],       46, 1, 1); 
      book(_h["el_eta_0j_pb"],        51, 1, 1);  
      book(_h["el_eta_1j_pb"],        56, 1, 1);   

      book(_h["Wplus_N_incl_pb"],      3, 1, 1);
      book(_h["Wplus_HT_1j_fb"],       8, 1, 1);
      book(_h["Wplus_W_pt_1j_fb"],    13, 1, 1);
      book(_h["Wplus_jet_pt1_1j_fb"], 18, 1, 1); 
      book(_h["Wplus_jet_y1_1j_fb"],  23, 1, 1);
      book(_h["Wplus_HT_2j_fb"],      38, 1, 1);
      book(_h["Wplus_W_pt_2j_fb"],    43, 1, 1);
      book(_h["Wplus_jet_pt1_2j_fb"], 48, 1, 1); 
      book(_h["Wplus_el_eta_0j_pb"],  53, 1, 1);
      book(_h["Wplus_el_eta_1j_pb"],  58, 1, 1);

      book(_h["Wminus_N_incl_pb"],     3, 1, 2);
      book(_h["Wminus_HT_1j_fb"],      8, 1, 2);
      book(_h["Wminus_W_pt_1j_fb"],   13, 1, 2);
      book(_h["Wminus_jet_pt1_1j_fb"],18, 1, 2);  
      book(_h["Wminus_jet_y1_1j_fb"], 23, 1, 2);
      book(_h["Wminus_HT_2j_fb"],     38, 1, 2);
      book(_h["Wminus_W_pt_2j_fb"],   43, 1, 2);
      book(_h["Wminus_jet_pt1_2j_fb"],48, 1, 2); 
      book(_h["Wminus_el_eta_0j_pb"], 53, 1, 2);
      book(_h["Wminus_el_eta_1j_pb"], 58, 1, 2);

      // and ratios
      book(_r["WplusOverWminus_N_incl_pb"],      3, 1, 3);
      book(_r["WplusOverWminus_HT_1j_fb"],       8, 1, 3);
      book(_r["WplusOverWminus_W_pt_1j_fb"],    13, 1, 3);
      book(_r["WplusOverWminus_jet_pt1_1j_fb"], 18, 1, 3);
      book(_r["WplusOverWminus_jet_y1_1j_fb"],  23, 1, 3);
      book(_r["WplusOverWminus_HT_2j_fb"],      38, 1, 3);
      book(_r["WplusOverWminus_W_pt_2j_fb"],    43, 1, 3);
      book(_r["WplusOverWminus_jet_pt1_2j_fb"], 48, 1, 3);
      book(_r["WplusOverWminus_el_eta_0j_pb"],  53, 1, 3);
      book(_r["WplusOverWminus_el_eta_1j_pb"],  58, 1, 3);

    }


    // Perform the per-event analysis
    void analyze(const Event& event) {
    
      // retrieve the dressed electrons
      const Particles& signal_leptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      if (signal_leptons.size() != 1 ) vetoEvent;
      const Particle& lepton = signal_leptons[0];

      // calulate MET and mT
      const Particles& invisibles = apply<InvisibleFinalState>(event, "InvFS").particles();
      FourMomentum pMET = sum(invisibles, Kin::mom, FourMomentum()).setZ(0);
      const double MET = pMET.pT() / GeV;
      const double mT = sqrt( 2 * lepton.Et() / GeV * MET * (1 - cos(deltaPhi(lepton,pMET))));
      if ( MET <= 25. ) vetoEvent;
      if ( mT  <= 40. ) vetoEvent;

      // retrieve jets
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
   
      idiscardIfAnyDeltaRLess(jets, signal_leptons, 0.2);

      // apply event selection on dR
      for (const Jet& j : jets) {
          if (deltaR(j, lepton) < 0.4) vetoEvent;
      }

      // calculate the observables
      const double w_pt = (lepton.momentum() + pMET).pT() / GeV;   
      const size_t njets = jets.size();
      double ST = sum(jets, Kin::pT, 0.0); // scalar pT sum of all selected jets
      const double HT = ST + lepton.pT() / GeV + MET; //missET;

      // fill W histograms		
      _h["N_pb"]->fill(njets);
      for (size_t i = 0; i <= njets; ++i) {
          _h["N_incl_pb"]->fill(i);
      }
      _h["el_eta_0j_pb"]->fill(lepton.abseta());

      if (njets > 0) {
          _h["HT_1j_fb"]->fill(HT);
          _h["W_pt_1j_fb"]->fill(w_pt);
          _h["jet_pt1_1j_fb"]->fill(jets[0].pT()/GeV);
          _h["jet_y1_1j_fb"]->fill(jets[0].absrap());
          _h["el_eta_1j_pb"]->fill(lepton.abseta());
      }
      if (njets > 1) {
          _h["HT_2j_fb"]->fill(HT);
          _h["W_pt_2j_fb"]->fill(w_pt);
          _h["jet_pt1_2j_fb"]->fill(jets[0].pT()/GeV);
          _h["DeltaRj12_2j_fb"]->fill(deltaR(jets[0],jets[1]));
          _h["jet_pt2_2j_fb"]->fill(jets[1].pT()/GeV);
          _h["jet_y2_2j_fb"]->fill(jets[1].absrap());
          _h["jet_mass12_2j_fb"]->fill( (jets[0].mom()+jets[1].mom()).mass()/GeV);
      }		
        // fill W+ histograms
      if (lepton.charge() > 0) {
          for (size_t i = 0; i <= njets; ++i) {
              _h["Wplus_N_incl_pb"]->fill(i);
          }
          _h["Wplus_el_eta_0j_pb"]->fill(lepton.abseta());
          if (njets > 0) {
              _h["Wplus_HT_1j_fb"]->fill(HT);
              _h["Wplus_W_pt_1j_fb"]->fill(w_pt);
              _h["Wplus_jet_pt1_1j_fb"]->fill(jets[0].pT()/GeV);
              _h["Wplus_jet_y1_1j_fb"]->fill(jets[0].absrap());
              _h["Wplus_el_eta_1j_pb"]->fill(lepton.abseta());
              if (njets > 1) {
                  _h["Wplus_HT_2j_fb"]->fill(HT);
                  _h["Wplus_W_pt_2j_fb"]->fill(w_pt);
                  _h["Wplus_jet_pt1_2j_fb"]->fill(jets[0].pT()/GeV);
              }		
          }
      }
      // fill W- histograms
      if (lepton.charge() < 0) {
        for (size_t i = 0; i <= njets; ++i) {
          _h["Wminus_N_incl_pb"]->fill(i);
        }
        _h["Wminus_el_eta_0j_pb"]->fill(lepton.abseta());
        if (njets > 0) {
          _h["Wminus_HT_1j_fb"]->fill(HT);
          _h["Wminus_W_pt_1j_fb"]->fill(w_pt);
          _h["Wminus_jet_pt1_1j_fb"]->fill(jets[0].pT()/GeV);
          _h["Wminus_jet_y1_1j_fb"]->fill(jets[0].absrap());
          _h["Wminus_el_eta_1j_pb"]->fill(lepton.abseta());
          if (njets > 1) {
            _h["Wminus_HT_2j_fb"]->fill(HT);
            _h["Wminus_W_pt_2j_fb"]->fill(w_pt);
            _h["Wminus_jet_pt1_2j_fb"]->fill(jets[0].pT()/GeV);
          }		
        }
      }
    }
    


    void finalize() {
      const double scalefactor_fb = crossSection() / sumOfWeights() / femtobarn;
      const double scalefactor_pb = crossSection() / sumOfWeights() / picobarn;

      for (auto& hit : _h){
        if (hit.first.find("_fb") != string::npos) scale(hit.second, scalefactor_fb);
        else                                       scale(hit.second, scalefactor_pb);
      }
      for (auto& rit : _r) {
        string ratio_label = "WplusOverWminus";
        string num_name   = rit.first;
        string denom_name = rit.first;
        num_name.replace(rit.first.find(ratio_label),ratio_label.length(),"Wplus");
        denom_name.replace(rit.first.find(ratio_label),ratio_label.length(),"Wminus");                   
        divide(_h[num_name], _h[denom_name], rit.second);
      }
    }

    protected:
      // Data members like post-cuts event weight counters go here
      size_t _mode;

    private:
      map<string, Histo1DPtr> _h;
      map<string, Scatter2DPtr> _r;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1635273);
}
