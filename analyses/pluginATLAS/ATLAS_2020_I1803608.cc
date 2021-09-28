// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  
  /// VBFZ in pp at 13 TeV
  class ATLAS_2020_I1803608 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2020_I1803608);

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 5.0);

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muons(Cuts::abspid == PID::MUON);

      Cut cuts_el = (Cuts::pT > 25*GeV) && ( Cuts::abseta < 1.37 || (Cuts::abseta > 1.52 && Cuts::abseta < 2.47) );
      Cut cuts_mu = (Cuts::pT > 25*GeV) && (Cuts::abseta < 2.4);

      DressedLeptons dressed_electrons(photons, electrons, 0.1, cuts_el);
      declare(dressed_electrons, "DressedElectrons");

      DressedLeptons dressed_muons(photons, muons, 0.1, cuts_mu);
      declare(dressed_muons, "DressedMuons");

      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "Jets");

      _doControl = bool(getOption("TYPE") != "EW_ONLY");
      if (_doControl) {
        initialisePlots(SRplots, "SR");
        initialisePlots(CRAplots, "CRA");
        initialisePlots(CRBplots, "CRB");
        initialisePlots(CRCplots, "CRC");
      }
      else {
        initialisePlots(SRplots, "EW");
     }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // access fiducial electrons and muons
      const Particle *l1 = nullptr, *l2 = nullptr;
      auto muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons();
      auto elecs = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons();

      // Dilepton selection 1: =2 leptons of the same kind
      if (muons.size()+elecs.size()!=2) vetoEvent;
      if      (muons.size()==2) { l1=&muons[0]; l2=&muons[1]; }
      else if (elecs.size()==2) { l1=&elecs[0]; l2=&elecs[1]; }
      else vetoEvent;

      // Dilepton selection 2: oppostie-charge and in mass range
      if ( !oppCharge(*l1, *l2) )  vetoEvent;
      if ( !inRange((l1->mom()+l2->mom()).mass()/GeV, 81.0, 101.0) ) vetoEvent;

      // Electron-jet overlap removal (note: muons are not included in jet finding)
      // make sure jets do not overlap with an electron within DR<0.2
      Jets jets;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 4.4)) {
        if (elecs.size() == 2 && (deltaR(j, *l1, RAPIDITY) < 0.2 || deltaR(j, *l2, RAPIDITY) < 0.2 )) {
          continue;
        }
        jets += j;
      }

      // Require 2 jets with pT > 85 and 80 GeV
      if (jets.size() < 2) vetoEvent;

      // Calculate the observables
      Variables vars(jets, l1, l2);

      // make sure neither lepton overlaps with a jet within 0.4
      for (const Jet& j : jets) {
        if (deltaR(j, *l1, RAPIDITY) < 0.4 || deltaR(j, *l2, RAPIDITY) < 0.4)  vetoEvent;
      }

      if (jets[0].pt() < 85*GeV || jets[1].pt() < 80*GeV ) vetoEvent;

      bool pass_VBFZtopo = (vars.mjj > 250*GeV && vars.Dyjj > 2.0 && vars.pTll > 20*GeV && vars.Zcent < 1.0 && vars.pTbal < 0.15);

      if (pass_VBFZtopo) {
        if      (_doControl && vars.Ngj  > 0 && vars.Zcent <  0.5) fillPlots(vars, CRAplots);
        else if (_doControl && vars.Ngj  > 0 && vars.Zcent >= 0.5) fillPlots(vars, CRBplots);
        else if (_doControl && vars.Ngj == 0 && vars.Zcent >= 0.5) fillPlots(vars, CRCplots);
        else {
          fillPlots(vars, SRplots);
        }
      }
    }

    void finalize() {
      const double xsec = crossSectionPerEvent()/femtobarn;
      scalePlots(SRplots, xsec);
      scalePlots(CRAplots, xsec);
      scalePlots(CRBplots, xsec);
      scalePlots(CRCplots, xsec);
    }

    struct Variables {
      Variables(const vector<Jet>& jets, const Particle* l1, const Particle* l2) {
        // get the jets
        assert(jets.size()>=2);
        FourMomentum j1 = jets[0].mom(), j2 = jets[1].mom();
        pTj1 = j1.pT()/GeV; pTj2 = j2.pT()/GeV;
        assert(pTj1 >= pTj2);
        
        // build dilepton system
        FourMomentum ll = (l1->mom() + l2->mom());
        pTll = ll.pT(); mll = ll.mass();
        
        Nj = jets.size();
        Dyjj = std::abs(j1.rap() - j2.rap());
        mjj = (j1 + j2).mass();
        Dphijj = ( j1.rap() > j2.rap() ) ? mapAngleMPiToPi(j1.phi() - j2.phi()) : mapAngleMPiToPi(j2.phi() - j1.phi());
        
        Jets gjets = getGapJets(jets);
        Ngj = gjets.size();
        pTgj = Ngj? gjets[0].pT()/GeV : 0;
        
        FourMomentum vecSum = (j1 + j2 + l1->mom() + l2->mom());
        double HT = j1.pT() + j2.pT() + l1->pT() + l2->pT();
        if (Ngj) { 
          vecSum += gjets[0].mom(); 
          HT += pTgj;
        }
        pTbal = vecSum.pT() / HT;
        
        Zcent = std::abs(ll.rap() - (j1.rap() + j2.rap())/2) / Dyjj;
      }
      
      double Zcent, pTj1, pTj2, pTgj, pTll, mll, Dyjj, mjj, Dphijj, pTbal;
      size_t Nj, Ngj;

      Jets getGapJets(const Jets& jets) {
        Jets gjets;
        if (jets.size() <= 2)  return gjets;
        FourMomentum j1 = jets[0].mom(), j2 = jets[1].mom();
        double yFwd = j1.rap(), yBwd = j2.rap();
        if (yFwd < yBwd) std::swap(yFwd,yBwd);
        for (size_t i = 2; i < jets.size(); ++i)
         if (inRange(jets[i].rap(), yBwd, yFwd)) gjets += jets[i];
        return gjets;
      }

    }; // struct variables

    struct Plots {
      string label;
      Histo1DPtr m_jj, Dphi_jj, Dy_jj, pT_ll;
    };

    void initialisePlots(Plots& plots, const string& phase_space) {
      plots.label   = phase_space;
      size_t region = 0;
      if (phase_space == "SR")   region = 4;
      if (phase_space == "CRA")  region = 8;
      if (phase_space == "CRB")  region = 12;
      if (phase_space == "CRC")  region = 16;
      book(plots.m_jj,    1 + region, 1, 1);
      book(plots.Dy_jj,   2 + region, 1, 1);
      book(plots.pT_ll,   3 + region, 1, 1);
      book(plots.Dphi_jj, 4 + region, 1, 1);
    }

    void fillPlots(const Variables& vars, Plots& plots) {
      // The mjj plot extends down to 250 GeV
      plots.m_jj->fill(vars.mjj/GeV);
      if (vars.mjj > 1000*GeV) {
        plots.Dy_jj->fill(vars.Dyjj);
        plots.Dphi_jj->fill(vars.Dphijj);
        plots.pT_ll->fill(vars.pTll/GeV);
      }
    }

    void scalePlots(Plots& plots, const double xsec){
      scale(plots.m_jj, xsec);
      scale(plots.Dy_jj, xsec);
      scale(plots.Dphi_jj, xsec);
      scale(plots.pT_ll, xsec);
    }

    //@}
    private:

      Plots SRplots, CRAplots, CRBplots, CRCplots;
      bool _doControl;
  };
  
  DECLARE_RIVET_PLUGIN(ATLAS_2020_I1803608);
}
