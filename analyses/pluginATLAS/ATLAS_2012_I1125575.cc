// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// ATLAS charged particle jet underlying event and jet radius dependence
  class ATLAS_2012_I1125575 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2012_I1125575()
      : Analysis("ATLAS_2012_I1125575")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState jet_input((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  0.5*GeV));
      declare(jet_input, "JET_INPUT");

      const ChargedFinalState track_input((Cuts::etaIn(-1.5, 1.5) && Cuts::pT >=  0.5*GeV));
      declare(track_input, "TRACK_INPUT");

      const FastJets jets02(jet_input, FastJets::ANTIKT, 0.2);
      declare(jets02, "JETS_02");

      const FastJets jets04(jet_input, FastJets::ANTIKT, 0.4);
      declare(jets04, "JETS_04");

      const FastJets jets06(jet_input, FastJets::ANTIKT, 0.6);
      declare(jets06, "JETS_06");

      const FastJets jets08(jet_input, FastJets::ANTIKT, 0.8);
      declare(jets08, "JETS_08");

      const FastJets jets10(jet_input, FastJets::ANTIKT, 1.0);
      declare(jets10, "JETS_10");

      // Mean number of tracks
      initializeProfiles(_h_meanNch, 1);

      // Mean of the average track pT in each region
      initializeProfiles(_h_meanPtAvg, 2);

      // Mean of the scalar sum of track pT in each region
      initializeProfiles(_h_meanPtSum, 3);

      // Distribution of Nch, in bins of leading track-jet pT
      initializeHistograms(_h_Nch, 4);

      // Distribution of average track-jet pT, in bins of leading track-jet pT
      initializeHistograms(_h_PtAvg, 5);

      // Distribution of sum of track-jet pT, in bins of leading track-jet pT
      initializeHistograms(_h_PtSum, 6);

      for (int i = 0; i < 5; ++i)
        book(_nEvents[i], "nEvents_"+to_str(i));
    }


    void initializeProfiles(Profile1DPtr plots[5][2], int distribution) {
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 2; ++j) {
          book(plots[i][j] ,distribution, i+1, j+1);
        }
      }
    }


    void initializeHistograms(BinnedHistogram plots[5][2], int distribution) {
      Scatter2D refscatter = refData(1, 1, 1);
      for (int i = 0; i < 5; ++i) {
        for (int y = 0; y < 2; ++y) {
          for (size_t j = 0; j < refscatter.numPoints(); ++j) {
            int histogram_number = ((j+1)*2)-((y+1)%2);
            double low_edge = refscatter.point(j).xMin();
            double high_edge = refscatter.point(j).xMax();
            Histo1DPtr tmp;
            plots[i][y].add(low_edge, high_edge, book(tmp, distribution, i+1, histogram_number));
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      vector<Jets*> all_jets;
      Jets jets_02 = apply<FastJets>(event, "JETS_02").jetsByPt(Cuts::pT > 4*GeV && Cuts::abseta < 1.5);
      all_jets.push_back(&jets_02);
      Jets jets_04 = apply<FastJets>(event, "JETS_04").jetsByPt(Cuts::pT > 4*GeV && Cuts::abseta < 1.5);
      all_jets.push_back(&jets_04);
      Jets jets_06 = apply<FastJets>(event, "JETS_06").jetsByPt(Cuts::pT > 4*GeV && Cuts::abseta < 1.5);
      all_jets.push_back(&jets_06);
      Jets jets_08 = apply<FastJets>(event, "JETS_08").jetsByPt(Cuts::pT > 4*GeV && Cuts::abseta < 1.5);
      all_jets.push_back(&jets_08);
      Jets jets_10 = apply<FastJets>(event, "JETS_10").jetsByPt(Cuts::pT > 4*GeV && Cuts::abseta < 1.5);
      all_jets.push_back(&jets_10);

      // Count the number of tracks in the away and transverse regions, for each set of jets
      double n_ch[5][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
      // Also add up the sum pT
      double sumpt[5][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
      // ptmean = sumpt / n_ch
      double ptavg[5][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
      // lead jet pT defines which bin we want to fill
      double lead_jet_pts[5] = {0.0};

      // Loop over each of the jet radii:
      for (int i = 0; i < 5; ++i) {
        if (all_jets[i]->size() < 1) continue;

        // Find the lead jet pT
        lead_jet_pts[i] = all_jets[i]->at(0).pT();

        // Loop over each of the charged particles
        const Particles& tracks = apply<ChargedFinalState>(event, "TRACK_INPUT").particlesByPt();
        for(const Particle& t : tracks) {

          // Get the delta-phi between the track and the leading jet
          double dphi = deltaPhi(all_jets[i]->at(0), t);

          // Find out which region this puts it in.
          // 0 = away region, 1 = transverse region, 2 = toward region
          int region = region_index(dphi);

          // If the track is in the toward region, ignore it.
          if (region == 2) continue;

          // Otherwise, increment the relevant counters
          ++n_ch[i][region];
          sumpt[i][region] += t.pT();

        }
        // Calculate the pT_avg for the away and transverse regions.
        // (And make sure we don't try to divide by zero.)
        ptavg[i][0] = (n_ch[i][0] == 0 ? 0.0 : sumpt[i][0] / n_ch[i][0]);
        ptavg[i][1] = (n_ch[i][1] == 0 ? 0.0 : sumpt[i][1] / n_ch[i][1]);

        _nEvents[i]->fill();
      }

      fillProfiles(_h_meanNch,    n_ch, lead_jet_pts, 1.0 / (2*PI));
      fillProfiles(_h_meanPtAvg, ptavg, lead_jet_pts, 1.0);
      fillProfiles(_h_meanPtSum, sumpt, lead_jet_pts, 1.0 / (2*PI));

      fillHistograms(_h_Nch,    n_ch, lead_jet_pts);
      fillHistograms(_h_PtAvg, ptavg, lead_jet_pts);
      fillHistograms(_h_PtSum, sumpt, lead_jet_pts);
    }


    void fillProfiles(Profile1DPtr plots[5][2], double var[5][2], double lead_pt[5], double scale) {
      for (int i=0; i<5; ++i) {
        double pt = lead_pt[i];
        for (int j=0; j<2; ++j) {
          double v = var[i][j];
          plots[i][j]->fill(pt, v*scale);
        }
      }
    }


    void fillHistograms(BinnedHistogram plots[5][2], double var[5][2], double lead_pt[5]) {
      for (int i=0; i<5; ++i) {
        double pt = lead_pt[i];
        for (int j=0; j<2; ++j) {
          double v = var[i][j];
          plots[i][j].fill(pt, v);
        }
      }
    }


    int region_index(double dphi) {
      assert(inRange(dphi, 0.0, PI, CLOSED, CLOSED));
      if (dphi < PI/3.0) return 2;
      if (dphi < 2*PI/3.0) return 1;
      return 0;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      finalizeHistograms(_h_Nch);
      finalizeHistograms(_h_PtAvg);
      finalizeHistograms(_h_PtSum);
    }


    void finalizeHistograms(BinnedHistogram plots[5][2]) {
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 2; ++j) {
          vector<Histo1DPtr> histos = plots[i][j].histos();
          for(Histo1DPtr h : histos) {
            scale(h, 1.0/ *_nEvents[i]);
          }
        }
      }
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    CounterPtr _nEvents[5];

    Profile1DPtr _h_meanNch[5][2];
    Profile1DPtr _h_meanPtAvg[5][2];
    Profile1DPtr _h_meanPtSum[5][2];

    BinnedHistogram _h_Nch[5][2];
    BinnedHistogram _h_PtAvg[5][2];
    BinnedHistogram _h_PtSum[5][2];

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1125575);

}
