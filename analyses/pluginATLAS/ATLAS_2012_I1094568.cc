// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {

  /// Top pair production with central jet veto
  class ATLAS_2012_I1094568 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2012_I1094568);

    struct Plots {
      // Track which veto region this is, to match the autobooked histograms
      int region_index;

      // Lower rapidity boundary or veto region
      double y_low;
      // Upper rapidity boundary or veto region
      double y_high;

      double veto_Q0;
      double veto_Qsum;

      // Histograms to store the veto jet pT and sum(veto jet pT) histograms.
      Histo1DPtr h_veto_Q0;
      Histo1DPtr h_veto_Qsum;

      // Scatter2Ds for the gap fractions
      Scatter2DPtr gapFrac_Q0;
      Scatter2DPtr gapFrac_Qsum;
    };

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 4.5);

      /// Get electrons from truth record
      FinalState elec_fs(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      FinalState muon_fs(Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      declare(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 4.5);
      neutrino_fs.acceptNeutrinos();
      declare(neutrino_fs, "NEUTRINO_FS");

      // Get the jets
      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(fs, "jet_input");
      declare(jets, "JETS");

      // get b-hadrons
      declare(HeavyHadrons(Cuts::pT > 5*GeV), "BHadrons");

      // Initialise weight counter
      book(m_total_weight, "_total_weight");

      // Init histogramming for the various regions
      m_plots[0].region_index = 1;
      m_plots[0].y_low = 0.0;
      m_plots[0].y_high = 0.8;
      initializePlots(m_plots[0]);
      //
      m_plots[1].region_index = 2;
      m_plots[1].y_low = 0.8;
      m_plots[1].y_high = 1.5;
      initializePlots(m_plots[1]);
      //
      m_plots[2].region_index = 3;
      m_plots[2].y_low = 1.5;
      m_plots[2].y_high = 2.1;
      initializePlots(m_plots[2]);
      //
      m_plots[3].region_index = 4;
      m_plots[3].y_low = 0.0;
      m_plots[3].y_high = 2.1;
      initializePlots(m_plots[3]);
    }


    void initializePlots(Plots& plots) {
      plots.veto_Q0 = 0.0;
      const string veto_Q0_name = "TMP/vetoJetPt_Q0_" + to_str(plots.region_index);
      book(plots.h_veto_Q0, veto_Q0_name, 200, 0.0, 1000.0);
      book(plots.gapFrac_Q0, plots.region_index, 1, 1, true);

      plots.veto_Qsum = 0.0;
      const string veto_Qsum_name = "TMP/vetoJetPt_Qsum_" + to_str(plots.region_index);
      book(plots.h_veto_Qsum, veto_Qsum_name, 200, 0.0, 1000.0);
      book(plots.gapFrac_Qsum, plots.region_index, 2, 1, true);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// Get the various sets of final state particles
      const Particles& elecFS = apply<FinalState>(event, "ELEC_FS").particlesByPt();
      const Particles& muonFS = apply<FinalState>(event, "MUON_FS").particlesByPt();
      const Particles& neutrinoFS = apply<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();

      // Get all jets with pT > 25 GeV and |y| < 2.4
      Jets jets = apply<FastJets>(event, "JETS").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 2.4);

      // For each of the jets that pass the rapidity cut, only keep those that are not
      // too close to any leptons
      idiscardIfAnyDeltaRLess(jets, elecFS, 0.4);
      idiscardIfAnyDeltaRLess(jets, muonFS, 0.4);

      // Get b hadrons with pT > 5 GeV
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();

      // For each of the good jets, check whether any are b-jets (via dR matching)
      size_t nMatches = 0;
      Jets bJets, vetoJets;
      for (const Jet& jet : jets) {
        bool isBjet = any(bHadrons, DeltaRLess(jet, 0.3));
        if (isBjet) { ++nMatches; bJets += jet; }
        if (!isBjet || nMatches > 2)  vetoJets += jet;
      }

      // Get the MET by taking the vector sum of all neutrinos
      /// @todo Use MissingMomentum instead?
      double MET = 0;
      FourMomentum p_MET;
      for(const Particle& p: neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();

      // Now we have everything we need to start doing the event selections
      bool passed_ee = false;

      // We want exactly 2 electrons...
      if (elecFS.size() == 2) {
        // ... with opposite sign charges.
        if (charge(elecFS[0]) != charge(elecFS[1])) {
          // Check the MET
          if (MET >= 40*GeV) {
            // Do some dilepton mass cuts
            const double dilepton_mass = (elecFS[0].momentum() + elecFS[1].momentum()).mass();
            if (dilepton_mass >= 15*GeV) {
              if (fabs(dilepton_mass - 91.0*GeV) >= 10.0*GeV) {
                // We need at least 2 b-jets
                passed_ee = bJets.size() > 1;
              }
            }
          }
        }
      }

      bool passed_mumu = false;
      // Now do the same checks for the mumu channel
      // So we now want 2 good muons...
      if (muonFS.size() == 2) {
        // ...with opposite sign charges.
        if (charge(muonFS[0]) != charge(muonFS[1])) {
          // Check the MET
          if (MET >= 40*GeV) {
            // and do some di-muon mass cuts
            const double dilepton_mass = (muonFS.at(0).momentum() + muonFS.at(1).momentum()).mass();
            if (dilepton_mass >= 15*GeV) {
              if (fabs(dilepton_mass - 91.0*GeV) >= 10.0*GeV) {
                // Need at least 2 b-jets
                passed_mumu = bJets.size() > 1;
              }
            }
          }
        }
      }

      bool passed_emu = false;
      // Finally, the same again with the emu channel
      // We want exactly 1 electron and 1 muon
      if (elecFS.size() == 1 && muonFS.size() == 1) {
        // With opposite sign charges
        if (charge(elecFS[0]) != charge(muonFS[0])) {
          // Calculate HT: scalar sum of the pTs of the leptons and all good jets
          double HT = sum(jets, pT, 0.);
          HT += elecFS[0].pT();
          HT += muonFS[0].pT();
          // Keep events with HT > 130 GeV
          if (HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            passed_emu = bJets.size() > 1;
          }
        }
      }

      if (passed_ee || passed_mumu || passed_emu) {
        // If the event passes the selection, we use it for all gap fractions
        m_total_weight->fill();

        // Loop over each veto jet
        for (const Jet& j : vetoJets) {
          const double pt = j.pT();
          const double rapidity = j.absrap();
          // Loop over each region
          for (size_t i = 0; i < 4; ++i) {
            // If the jet falls into this region, get its pT and increment sum(pT)
            if (inRange(rapidity, m_plots[i].y_low, m_plots[i].y_high)) {
              m_plots[i].veto_Qsum += pt;
              // If we've already got a veto jet, don't replace it
              if (m_plots[i].veto_Q0 == 0.0)  m_plots[i].veto_Q0 = pt;
            }
          }
        }
        for (size_t i = 0; i < 4; ++i) {
          m_plots[i].h_veto_Q0->fill(m_plots[i].veto_Q0);
          m_plots[i].h_veto_Qsum->fill(m_plots[i].veto_Qsum);
          m_plots[i].veto_Q0 = 0.0;
          m_plots[i].veto_Qsum = 0.0;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double totalWeight = m_total_weight->val();
      for (size_t i = 0; i < 4; ++i) {
        finalizeGapFraction(totalWeight, m_plots[i].gapFrac_Q0,   m_plots[i].h_veto_Q0);
        finalizeGapFraction(totalWeight, m_plots[i].gapFrac_Qsum, m_plots[i].h_veto_Qsum);
      }
    }


    /// Convert temporary histos to cumulative efficiency scatters
    /// @todo Should be possible to replace this with a couple of YODA one-lines for diff -> integral and "efficiency division"
    void finalizeGapFraction(const double total_weight, Scatter2DPtr gapFrac, Histo1DPtr vetoPt) {
      // Stores the cumulative frequency of the veto jet pT histogram
      double vetoPtWeightSum = 0.0;

      // Keep track of which gap fraction point we're currently populating (#final_points != #tmp_bins)
      size_t fgap_point = 0;
      for (size_t i = 0; i < vetoPt->numBins(); ++i) {
        // If we've done the last "final" point, stop
        if (fgap_point == gapFrac->numPoints())  break;

        // Increment the cumulative vetoPt counter for this temp histo bin
        /// @todo Get rid of this and use vetoPt->integral(i+1) when points and bins line up?
        vetoPtWeightSum += vetoPt->bin(i).sumW();

        // If this temp histo bin's upper edge doesn't correspond to the reference point, don't finalise the scatter.
        // Note that points are ON the bin edges and have no width: they represent the integral up to exactly that point.
        if ( !fuzzyEquals(vetoPt->bin(i).xMax(), gapFrac->point(fgap_point).x()) )  continue;

        // Calculate the gap fraction and its uncertainty
        const double frac = (total_weight != 0.0) ? vetoPtWeightSum/total_weight : 0;
        const double fracErr = (total_weight != 0.0) ? sqrt(frac*(1-frac)/total_weight) : 0;
        gapFrac->point(fgap_point).setY(frac, fracErr);

        ++fgap_point;
      }
    }


  private:

    CounterPtr m_total_weight;
    Plots m_plots[4];

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094568);
}
