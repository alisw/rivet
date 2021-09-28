// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"


namespace Rivet {


  /// @brief WW production at 13 TeV
  class ATLAS_2021_I1852328 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2021_I1852328);


    /// Book histograms and initialise projections before the run
    void init() {
     
      const FinalState fs(Cuts::abseta < 5.);

      // Project photons for dressing
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Cuts for leptons
      Cut lepton_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 27*GeV);
          
      // Project dressed leptons (e/mu not from tau) with pT > 27 GeV and |eta| < 2.5
      PromptFinalState lep_bare(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      lep_bare.acceptTauDecays(false);
      DressedLeptons lep_dressed(photons, lep_bare, 0.1, lepton_cuts, false);
      declare(lep_dressed,"lep_dressed");

      // Leptons (+ dressed photons) to be removed from jets
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true); // true = use muons from prompt tau decays
      DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, Cuts::abseta < 2.5, true);
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true); // true = use electrons from prompt tau decays
      DressedLeptons all_dressed_el(photons, bare_el, 0.1, Cuts::abseta < 2.5, true);

      // Define hadrons as everything but dressed leptons (for jet clustering)
      VetoedFinalState hadrons(fs);
      hadrons.addVetoOnThisFinalState(all_dressed_el);
      hadrons.addVetoOnThisFinalState(all_dressed_mu);
      declare(hadrons, "hadrons");

      // Get MET from generic invisibles
      MissingMomentum mm(fs);
      declare(mm, "met");
      
      // Project jets
      FastJets jets(hadrons, FastJets::ANTIKT, 0.4, FastJets::Muons::ALL, FastJets::Invisibles::DECAY);
      declare(jets, "jets");
    
      book(_h["xs_inf"], 1, 1, 1);
      book(_h["xs_bveto_inf"], 1, 1, 2);
      book(_h["lep0pt"], 2, 1, 1);
      book(_h["lep0pt_bveto"], 2, 1, 2);
      book(_h["lep1pt"], 5, 1, 1);
      book(_h["lep1pt_bveto"], 5, 1, 2);
      book(_h["jet0pt"], 8, 1, 1);
      book(_h["jet0pt_bveto"], 8, 1, 2);
      book(_h["htjet"], 11, 1, 1);
      book(_h["htjet_bveto"], 11, 1, 2);
      book(_h["st"], 14, 1, 1);
      book(_h["st_bveto"], 14, 1, 2);
      book(_h["mt"], 17, 1, 1);
      book(_h["mt_bveto"], 17, 1, 2);
      book(_h["mll_inf"], 20, 1, 1);
      book(_h["mll_bveto_inf"], 20, 1, 2);
      book(_h["ptll"], 23, 1, 1);
      book(_h["ptll_bveto"], 23, 1, 2);
      book(_h["dphill_inf"], 26, 1, 1);
      book(_h["dphill_bveto_inf"], 26, 1, 2);
      book(_h["yll_inf"], 29, 1, 1);
      book(_h["yll_bveto_inf"], 29, 1, 2);
      book(_h["costhetastar_inf"], 32, 1, 1);
      book(_h["costhetastar_bveto_inf"], 32, 1, 2);
      book(_h["njet"], 35, 1, 1);
      book(_h["njet_bveto"], 35, 1, 2);

      book(_h["mll_jet0pt200_inf"], 38, 1, 1);
      book(_h["mll_jet0pt200_bveto_inf"], 38, 1, 2);
      book(_h["dphill_jet0pt200_inf"], 41, 1, 1);
      book(_h["dphill_jet0pt200_bveto_inf"], 41, 1, 2);

      book(_h["dphil1j0_lep0pt200_inf"], 44, 1, 1);
      book(_h["dphil1j0_lep0pt200_bveto_inf"], 44, 1, 2);
      book(_h["drl1j0_lep0pt200_inf"], 47, 1, 1);
      book(_h["drl1j0_lep0pt200_bveto_inf"], 47, 1, 2);
      book(_h["rl1l0_lep0pt200_inf"], 50, 1, 1);
      book(_h["rl1l0_lep0pt200_bveto_inf"], 50, 1, 2);
      book(_h["rl1j0_lep0pt200"], 53, 1, 1);
      book(_h["rl1j0_lep0pt200_bveto"], 53, 1, 2);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get invisible
      const MissingMomentum& met = apply<MissingMomentum>(event, "met");

      // Find leptons and sort by pT
      const Particles& leptons = apply<ParticleFinder>(event, "lep_dressed").particlesByPt();

      const Jets& jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 4.5 && Cuts::pT > 30*GeV);
      const Jets& btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
      
      // Define observables
 
      const FourMomentum dilep = leptons.size()>1 ? leptons[0].momentum() + leptons[1].momentum() : FourMomentum(0,0,0,0);
      const double lep0pt = leptons.size()>0 ? leptons[0].pT()/GeV : -1;
      const double lep1pt = leptons.size()>1 ? leptons[1].pT()/GeV : -1;
      const double mll = leptons.size()>1 ? dilep.mass()/GeV : -1;
      const double ptll = leptons.size()>1 ? dilep.pT()/GeV : -1;
      const double yll = leptons.size()>1 ? dilep.absrap(): -1; 
      const double dphill = leptons.size()>1 ? fabs(deltaPhi(leptons[0], leptons[1])) : -1.;
      const double costhetastar = leptons.size()>1 ? fabs(tanh((leptons[0].eta() - leptons[1].eta()) / 2)) : -1;
      const double jet0pt = jets.size()>0 ? jets[0].pT()/GeV : -1;
      const double dphil1j0 = leptons.size()>1 and jets.size()>0 ? fabs(deltaPhi(leptons[1], jets[0])) : -1.;
      const double drl1j0 = leptons.size()>1 and jets.size()>0 ? fabs(deltaR(leptons[1], jets[0])) : -1.;
      const double rl1l0 = leptons.size()>1 ? leptons[1].pT()/leptons[0].pT() : -1;
      const double rl1j0 = leptons.size()>1 and jets.size()>0 ? leptons[1].pT()/jets[0].pT() : -1;
      Vector3 dilep_tvec(dilep.px(), dilep.py(),0);
      Vector3 met_tvec = met.vectorMissingPt();
      double et_ll = std::sqrt(dilep_tvec.mod2() + dilep.mass()*dilep.mass());
      const double mt =sqrt(std::pow(et_ll + met.met() ,2) - (dilep_tvec + met_tvec).mod2());
      const int njet = jets.size();
   
      size_t nbtags = count(btagging_jets, hasBTag());
      const double htjet = sum(jets, pT, 0.0);
      const double st = sum(leptons, pT, htjet);

      if (leptons.size()!=2 ) vetoEvent;
      if (leptons[0].abspid() == leptons[1].abspid()) vetoEvent;
      if (leptons[0].pid()*leptons[1].pid()>0) vetoEvent;
      if (dilep.mass() <= 85*GeV)  vetoEvent; 
      if (jets.empty()) vetoEvent;

      _h["xs_inf"]->fill(0.5);
      _h["lep0pt"]->fill(lep0pt);
      _h["lep1pt"]->fill(lep1pt);
      _h["mll_inf"]->fill(mll);
      _h["ptll"]->fill(ptll);
      _h["yll_inf"]->fill(yll);
      _h["dphill_inf"]->fill(dphill);
      _h["costhetastar_inf"]->fill(costhetastar);
      _h["jet0pt"]->fill(jet0pt);
      _h["mt"]->fill(mt);
      _h["njet"]->fill(njet);
      _h["htjet"]->fill(htjet/GeV);
      _h["st"]->fill(st/GeV);
      
      if (lep0pt>200){
	      _h["dphil1j0_lep0pt200_inf"]->fill(dphil1j0);
	      _h["drl1j0_lep0pt200_inf"]->fill(drl1j0);
	      _h["rl1l0_lep0pt200_inf"]->fill(rl1l0);
	      _h["rl1j0_lep0pt200"]->fill(rl1j0);
      }
      if(jet0pt>200){
	     _h["mll_jet0pt200_inf"]->fill(mll);
	     _h["dphill_jet0pt200_inf"]->fill(dphill);
      }
      

      if(nbtags == 0) {
        _h["xs_bveto_inf"]->fill(0.5);
        _h["lep0pt_bveto"]->fill(lep0pt);
        _h["lep1pt_bveto"]->fill(lep1pt);
        _h["mll_bveto_inf"]->fill(mll);
        _h["ptll_bveto"]->fill(ptll);
        _h["yll_bveto_inf"]->fill(yll);
        _h["dphill_bveto_inf"]->fill(dphill);
        _h["costhetastar_bveto_inf"]->fill(costhetastar);
        _h["jet0pt_bveto"]->fill(jet0pt);
        _h["mt_bveto"]->fill(mt);
        _h["njet_bveto"]->fill(njet);
        _h["htjet_bveto"]->fill(htjet);
        _h["st_bveto"]->fill(st);
      
         if (lep0pt>200){
        	  _h["dphil1j0_lep0pt200_bveto_inf"]->fill(dphil1j0);
 	          _h["drl1j0_lep0pt200_bveto_inf"]->fill(drl1j0);
	          _h["rl1l0_lep0pt200_bveto_inf"]->fill(rl1l0);
    	      _h["rl1j0_lep0pt200_bveto"]->fill(rl1j0);
           }
          
         if(jet0pt>200){
	          _h["mll_jet0pt200_bveto_inf"]->fill(mll);
	          _h["dphill_jet0pt200_bveto_inf"]->fill(dphill);
          }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    
       const double sf = crossSectionPerEvent()/femtobarn;
      
       scale(_h, sf);
       for (auto& hist : _h) {
          if (hist.first.find("_inf") != string::npos) {
            hist.second->fillBin(hist.second->numBins()-1, hist.second->overflow().sumW());
            hist.second->overflow().reset();
          }
       }
    }

     //@}

  private:

    /// @name Histograms
    map<string, Histo1DPtr> _h;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2021_I1852328);
}
