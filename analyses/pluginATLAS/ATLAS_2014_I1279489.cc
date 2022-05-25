// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  struct Plots {
    string label;

    Histo1DPtr h_dy;
    Histo1DPtr h_mjj;
    Histo1DPtr h_njets;
    Histo1DPtr h_dphijj;
    Histo1DPtr h_ptbal;

    Histo1DPtr h_jetveto_mjj_veto;
    Histo1DPtr h_jetveto_mjj_inc;
    Histo1DPtr h_jetveto_dy_veto;
    Histo1DPtr h_jetveto_dy_inc;

    Histo1DPtr h_ptbaleff_mjj_veto;
    Histo1DPtr h_ptbaleff_mjj_inc;
    Histo1DPtr h_ptbaleff_dy_veto;
    Histo1DPtr h_ptbaleff_dy_inc;

    Scatter2DPtr s_jetveto_mjj;
    Scatter2DPtr s_jetveto_dy;

    Scatter2DPtr s_ptbaleff_mjj;
    Scatter2DPtr s_ptbaleff_dy;

    Profile1DPtr p_avgnjets_dy;
    Profile1DPtr p_avgnjets_mjj;
  };


  struct Variables {

    Variables(const Jets& jets, const Particle* lep1, const Particle* lep2) {
      FourMomentum j1 = jets[0].momentum();
      FourMomentum j2 = jets[1].momentum();
      jet1pt = j1.pT();
      jet2pt = j2.pT();
      assert(jet1pt > jet2pt);

      zpt = (lep1->mom() + lep2->mom()).pT();

      deltay = fabs(j1.rapidity() - j2.rapidity());
      mjj = (j1 + j2).mass();
      deltaphijj = deltaPhi(j1, j2) / PI;

      FourMomentum gapjet(0., 0., 0., 0.);
      ngapjets = _getNumGapJets(jets, gapjet);

      double ptbal_vec = (j1 + j2 + lep1->mom() + lep2->mom()).pT();
      double ptbal_sc = j1.pT() + j2.pT() + lep1->pT() + lep2->pT();
      ptbalance2 = ptbal_vec / ptbal_sc;

      double ptbal3_vec = (j1 + j2 + gapjet + lep1->mom() + lep2->mom()).pT();
      double ptbal3_sc = j1.pT() + j2.pT() + gapjet.pT() + lep1->pT() + lep2->pT();
      ptbalance3 = ptbal3_vec / ptbal3_sc;

      pass_jetveto = gapjet.pT() < 25.0*GeV;
      pass_ptbaleff = ptbalance2 < 0.15;
    }


    double jet1pt;
    double jet2pt;
    double zpt;

    double deltay;
    double mjj;
    double deltaphijj;
    double ptbalance2;
    double ptbalance3;
    int ngapjets;

    double dilepton_dr;

    bool pass_jetveto;
    bool pass_ptbaleff;


  private:

    bool _isBetween(const Jet& probe, const Jet& boundary1, const Jet& boundary2) {
      double y_p = probe.rap();
      double y_b1 = boundary1.rap();
      double y_b2 = boundary2.rap();

      double y_min = std::min(y_b1, y_b2);
      double y_max = std::max(y_b1, y_b2);

      if (y_p > y_min && y_p < y_max) return true;
      else return false;
    }

    int _getNumGapJets(const Jets& jets, FourMomentum& thirdJet) {
      if (jets.size() < 2) return 0;
      int n_between = 0;
      // Start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); ++i) {
        // If this jet is between the boundary jets and is hard enough, increment counter
        if (_isBetween(jets[i], jets[0], jets[1])) {
          if (n_between == 0) thirdJet = jets[i].momentum();
          ++n_between;
        }
      }
      return n_between;
    }

  };



  class ATLAS_2014_I1279489 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1279489);


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(Cuts::abseta < 5.);

      FinalState photon_fs(fs, Cuts::abspid == PID::PHOTON);

      FinalState electron_fs(fs, Cuts::abspid == PID::ELECTRON);

      FinalState muon_fs(fs, Cuts::abspid == PID::MUON);

      DressedLeptons dressed_electrons(photon_fs, electron_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_electrons, "DressedElectrons");

      DressedLeptons dressed_muons(photon_fs, muon_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_muons, "DressedMuons");

      FastJets jets(fs, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      initialisePlots(baseline_plots, "baseline");
      initialisePlots(highpt_plots,   "highpt");
      initialisePlots(search_plots,   "search");
      initialisePlots(control_plots,  "control");
      initialisePlots(highmass_plots, "highmass");
    }


    void initialisePlots(Plots& plots, const string& phase_space){
      /****************************************
       * Plot labeling:                       *
       * format = d0_-x0_-y0_                 *
       * d01 = baseline fiducial region       *
       * d02 = high-pt fiducial region        *
       * d03 = search fiducial region         *
       * d04 = control fiducial region        *
       * d05 = high-mass fiducial region      *
       *                                      *
       * x01 = mjj on x-axis                  *
       * x02 = delta-y on x-axis              *
       * x03 = njets on x-axis                *
       * x04 = dphijj on x-axis               *
       * x05 = ptbalance on x-axis            *
       *                                      *
       * y01 = differential cross-section     *
       * y02 = jet veto efficiency            *
       * y03 = ptbalance efficiency           *
       * y04 = average njets                  *
       ****************************************/
      plots.label = phase_space;

      if (phase_space=="baseline") {
        book(plots.h_mjj ,1, 1, 1);
        book(plots.h_dy ,3, 1, 1);

        book(plots.h_jetveto_mjj_veto,"_jetveto_mjj_baseline_veto", refData(8,1,1));
        book(plots.h_jetveto_mjj_inc, "_jetveto_mjj_baseline_inc",  refData(8,1,1));
        book(plots.h_jetveto_dy_veto, "_jetveto_dy_baseline_veto",  refData(9,1,1));
        book(plots.h_jetveto_dy_inc,  "_jetveto_dy_baseline_inc",   refData(9,1,1));

        book(plots.h_ptbaleff_mjj_veto, "_ptbaleff_mjj_baseline_veto", refData(12,1,1));
        book(plots.h_ptbaleff_mjj_inc,  "_ptbaleff_mjj_baseline_inc",  refData(12,1,1));
        book(plots.h_ptbaleff_dy_veto,  "_ptbaleff_dy_baseline_veto",  refData(13,1,1));
        book(plots.h_ptbaleff_dy_inc,   "_ptbaleff_dy_baseline_inc",   refData(13,1,1));

        book(plots.s_jetveto_mjj,   8, 1, 1);
        book(plots.s_jetveto_dy,    9, 1, 1);
        book(plots.s_ptbaleff_mjj, 12, 1, 1);
        book(plots.s_ptbaleff_dy,  13, 1, 1);

        book(plots.p_avgnjets_mjj ,10,1,1);
        book(plots.p_avgnjets_dy ,11,1,1);
      }

      if (phase_space=="highpt") {
        book(plots.h_mjj , 14, 1, 1);
        book(plots.h_dy , 16, 1, 1);

        book(plots.h_jetveto_mjj_veto, "_jetveto_mjj_highpt_veto", refData(18,1,1));
        book(plots.h_jetveto_mjj_inc,  "_jetveto_mjj_highpt_inc",  refData(18,1,1));
        book(plots.h_jetveto_dy_veto,  "_jetveto_dy_highpt_veto",  refData(19,1,1));
        book(plots.h_jetveto_dy_inc,   "_jetveto_dy_highpt_inc",   refData(19,1,1));

        book(plots.h_ptbaleff_mjj_veto, "_ptbaleff_mjj_highpt_veto", refData(22,1,1));
        book(plots.h_ptbaleff_mjj_inc,  "_ptbaleff_mjj_highpt_inc",  refData(22,1,1));
        book(plots.h_ptbaleff_dy_veto,  "_ptbaleff_dy_highpt_veto",  refData(23,1,1));
        book(plots.h_ptbaleff_dy_inc,   "_ptbaleff_dy_highpt_inc",   refData(23,1,1));

        book(plots.s_jetveto_mjj,  18, 1, 1);
        book(plots.s_jetveto_dy,   19, 1, 1);
        book(plots.s_ptbaleff_mjj, 22, 1, 1);
        book(plots.s_ptbaleff_dy,  23, 1, 1);

        book(plots.p_avgnjets_mjj , 20,1,1);
        book(plots.p_avgnjets_dy , 21,1,1);
      }

      if (phase_space=="search") {
        book(plots.h_mjj , 2,1,1);
        book(plots.h_dy , 4,1,1);
      }

      if (phase_space=="control") {
        book(plots.h_mjj , 15,1,1);
        book(plots.h_dy , 17,1,1);
      }

      if (phase_space=="highmass") {
        book(plots.h_njets , 5,1,1);
        book(plots.h_dphijj , 7,1,1);
        book(plots.h_ptbal , 6,1,1);
      }
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Make sure that we have a Z-candidate:
      const Particle *lep1 = nullptr, *lep2 = nullptr;
      //
      const vector<DressedLepton>& muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons();
      if (muons.size() == 2) {
        const FourMomentum dimuon = muons[0].mom() + muons[1].mom();
        if ( inRange(dimuon.mass()/GeV, 81.0, 101.0) && PID::charge3(muons[0].pid()) != PID::charge3(muons[1].pid()) ) {
          lep1 = &muons[0];
          lep2 = &muons[1];
        }
      }
      //
      const vector<DressedLepton>& electrons = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons();
      if (electrons.size() == 2) {
        const FourMomentum dielectron = electrons[0].mom() + electrons[1].mom();
        if ( inRange(dielectron.mass()/GeV, 81.0, 101.0) && PID::charge3(electrons[0].pid()) != PID::charge3(electrons[1].pid()) ) {
          if (lep1 && lep2) {
            MSG_INFO("Found Z candidates using both electrons and muons! Continuing with the muon-channel candidate");
          } else {
            lep1 = &electrons[0];
            lep2 = &electrons[1];
          }
        }
      }
      // If there's no Z-candidate, we won't use this event:
      if (!lep1 || !lep2) vetoEvent;


      // Do lepton-jet overlap removal:
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 4.4);
      idiscardIfAnyDeltaRLess(jets, muons, 0.3);
      idiscardIfAnyDeltaRLess(jets, electrons, 0.3);

      // If we don't have at least 2 good jets, we won't use this event.
      if (jets.size() < 2) vetoEvent;


      // Plotting, using variables and histo classes calculated by the Variables object constructor
      Variables vars(jets, lep1, lep2);
      bool pass_baseline = (vars.jet1pt > 55*GeV && vars.jet2pt > 45*GeV);
      bool pass_highpt   = (vars.jet1pt > 85*GeV && vars.jet2pt > 75*GeV);
      bool pass_highmass = (pass_baseline && vars.mjj > 1000*GeV);
      bool pass_search   = (pass_baseline && vars.zpt > 20*GeV && vars.ngapjets == 0 && vars.ptbalance2 < 0.15 && vars.mjj > 250*GeV);
      bool pass_control  = (pass_baseline && vars.zpt > 20*GeV && vars.ngapjets  > 0 && vars.ptbalance3 < 0.15 && vars.mjj > 250*GeV);
      //
      if (pass_baseline) fillPlots(vars, baseline_plots, "baseline");
      if (pass_highpt)   fillPlots(vars, highpt_plots,   "highpt");
      if (pass_highmass) fillPlots(vars, highmass_plots, "highmass");
      if (pass_search)   fillPlots(vars, search_plots,   "search");
      if (pass_control)  fillPlots(vars, control_plots,  "control");
    }


    void fillPlots(const Variables& vars, Plots& plots, string phase_space) {
      if (phase_space == "baseline" || phase_space == "highpt" || phase_space == "search" || phase_space == "control") {
        plots.h_dy->fill(vars.deltay);
        plots.h_mjj->fill(vars.mjj);
      }

      if (phase_space == "baseline" || phase_space == "highpt") {
        if (vars.pass_jetveto) {
          plots.h_jetveto_dy_veto->fill(vars.deltay);
          plots.h_jetveto_mjj_veto->fill(vars.mjj);
        }
        plots.h_jetveto_dy_inc->fill(vars.deltay);
        plots.h_jetveto_mjj_inc->fill(vars.mjj);

        if (vars.pass_ptbaleff) {
          plots.h_ptbaleff_mjj_veto->fill(vars.mjj);
          plots.h_ptbaleff_dy_veto->fill(vars.deltay);
        }
        plots.h_ptbaleff_mjj_inc->fill(vars.mjj);
        plots.h_ptbaleff_dy_inc->fill(vars.deltay);

        plots.p_avgnjets_dy->fill(vars.deltay, vars.ngapjets);
        plots.p_avgnjets_mjj->fill(vars.mjj, vars.ngapjets);
      }

      if (phase_space == "highmass") {
        plots.h_njets->fill(vars.ngapjets);
        plots.h_dphijj->fill(vars.deltaphijj);
        plots.h_ptbal->fill(vars.ptbalance2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      finalizePlots(baseline_plots);
      finalizePlots(highpt_plots);
      finalizePlots(search_plots);
      finalizePlots(control_plots);
      finalizePlots(highmass_plots);
      finalizeEfficiencies(baseline_plots);
      finalizeEfficiencies(highpt_plots);
    }

    void finalizePlots(Plots& plots) {
      if (plots.h_dy)     normalize(plots.h_dy);
      if (plots.h_mjj)    normalize(plots.h_mjj);
      if (plots.h_dphijj) normalize(plots.h_dphijj);
      if (plots.h_njets)  normalize(plots.h_njets);
      if (plots.h_ptbal)  normalize(plots.h_ptbal);
    }

    void finalizeEfficiencies(Plots& plots) {

      if (plots.label != "baseline" && plots.label != "highpt")  return;
      //size_t offset = plots.label == "baseline"? 0 : 10;

      if (plots.h_jetveto_mjj_veto && plots.h_jetveto_mjj_inc) {
        divide(plots.h_jetveto_mjj_veto, plots.h_jetveto_mjj_inc, plots.s_jetveto_mjj);
      }
      //getScatter2D(8+offset, 1, 1)->addAnnotation("InclusiveSumWeights", plots.h_jetveto_mjj_inc->integral());

      if (plots.h_jetveto_dy_veto && plots.h_jetveto_dy_inc) {
        divide(plots.h_jetveto_dy_veto, plots.h_jetveto_dy_inc, plots.s_jetveto_dy);
      }

      if (plots.h_ptbaleff_mjj_veto && plots.h_ptbaleff_mjj_inc) {
        divide(plots.h_ptbaleff_mjj_veto, plots.h_ptbaleff_mjj_inc, plots.s_ptbaleff_mjj);
      }

      if (plots.h_ptbaleff_dy_veto && plots.h_ptbaleff_dy_inc) {
        divide(plots.h_ptbaleff_dy_veto, plots.h_ptbaleff_dy_inc, plots.s_ptbaleff_dy);
      }
    }

    //@}


  private:

    //Variables* vars;

    Plots baseline_plots;
    Plots highpt_plots;
    Plots search_plots;
    Plots control_plots;
    Plots highmass_plots;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1279489);

}
