// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

using namespace std;

namespace Rivet {

  // UE charged particles vs. leading jet
  class CMS_2015_PAS_FSQ_15_007 : public Analysis {
  public:
    /// Constructor
    CMS_2015_PAS_FSQ_15_007() : Analysis("CMS_2015_PAS_FSQ_15_007") {}

    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 2.0 && Cuts::pt > 500 * MeV);
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet(Cuts::abseta < 2.5 && Cuts::pt > 500 * MeV);
      const FastJets jetpro(cfsforjet, FastJets::SISCONE, 0.5);
      declare(jetpro, "Jets");

      book(_h_PtSum_vs_leadTrackPt_transMin, 1, 1, 1);
      book(_h_PtSum_vs_leadTrackPt_transMax, 2, 1, 1);
      book(_h_PtSum_vs_leadTrackPt_transDiff, 3, 1, 1);
      book(_h_PtSum_vs_leadTrackPt_transAvg, 4, 1, 1);

      book(_h_Nch_vs_leadTrackPt_transMin, 5, 1, 1);
      book(_h_Nch_vs_leadTrackPt_transMax, 6, 1, 1);
      book(_h_Nch_vs_leadTrackPt_transDiff, 7, 1, 1);
      book(_h_Nch_vs_leadTrackPt_transAvg, 8, 1, 1);

      book(_h_PtSum_vs_leadJetPt_transMin, 9, 1, 1);
      book(_h_PtSum_vs_leadJetPt_transMax, 10, 1, 1);
      book(_h_PtSum_vs_leadJetPt_transDiff, 11, 1, 1);
      book(_h_PtSum_vs_leadJetPt_transAvg, 12, 1, 1);

      book(_h_Nch_vs_leadJetPt_transMin, 13, 1, 1);
      book(_h_Nch_vs_leadJetPt_transMax, 14, 1, 1);
      book(_h_Nch_vs_leadJetPt_transDiff, 15, 1, 1);
      book(_h_Nch_vs_leadJetPt_transAvg, 16, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_leadjet, p_leadtrack;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(1.0 * GeV)) {
        if (j.abseta() < 2.0) {
          p_leadjet = j.momentum();
          break;
        }
      }

      for (const Particle& j : apply<ChargedFinalState>(event, "CFS").particlesByPt(0.5 * GeV)) {
        if (j.abseta() < 2.0) {
          p_leadtrack = j.momentum();
          break;
        }
      }

      if (p_leadjet.isZero() && p_leadtrack.isZero())
        vetoEvent;
      const double phileadjet = p_leadjet.phi();
      const double pTleadjet = p_leadjet.pT();

      const double phileadtrack = p_leadtrack.phi();
      const double pTleadtrack = p_leadtrack.pT();

      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();

      int nTransverse_leadjet = 0;
      double ptSumTransverse_leadjet = 0.;
      int nTransverse1_leadjet = 0;
      double ptSumTransverse1_leadjet = 0.;
      int nTransverse2_leadjet = 0;
      double ptSumTransverse2_leadjet = 0.;
      int nTransverseMin_leadjet = 0;
      double ptSumTransverseMin_leadjet = 0.;
      int nTransverseMax_leadjet = 0;
      double ptSumTransverseMax_leadjet = 0.;
      int nTowards_leadjet = 0;
      double ptSumTowards_leadjet = 0.;
      int nAway_leadjet = 0;
      double ptSumAway_leadjet = 0.;

      int nTransverse_leadtrack = 0;
      double ptSumTransverse_leadtrack = 0.;
      int nTransverse1_leadtrack = 0;
      double ptSumTransverse1_leadtrack = 0.;
      int nTransverse2_leadtrack = 0;
      double ptSumTransverse2_leadtrack = 0.;
      int nTransverseMin_leadtrack = 0;
      double ptSumTransverseMin_leadtrack = 0.;
      int nTransverseMax_leadtrack = 0;
      double ptSumTransverseMax_leadtrack = 0.;
      int nTowards_leadtrack = 0;
      double ptSumTowards_leadtrack = 0.;
      int nAway_leadtrack = 0;
      double ptSumAway_leadtrack = 0.;

      for (const Particle& p : particles) {
        const double pT = p.pT() / GeV;

        if (!p_leadjet.isZero()) {
          double dphi_leadjet = p.phi() - phileadjet;
          while (dphi_leadjet > PI) {
            dphi_leadjet = dphi_leadjet - 2.0 * PI;
          }
          while (dphi_leadjet < -PI) {
            dphi_leadjet = dphi_leadjet + 2. * PI;
          }

          if (dphi_leadjet > PI / 3. && dphi_leadjet < PI * 2. / 3.) {  // Transverse1 region
            nTransverse_leadjet++;
            ptSumTransverse_leadjet += pT;
            nTransverse1_leadjet++;
            ptSumTransverse1_leadjet += pT;
          }

          if (dphi_leadjet < -PI / 3. && dphi_leadjet > -PI * 2. / 3.) {  // Transverse2 region
            nTransverse_leadjet++;
            ptSumTransverse_leadjet += pT;
            nTransverse2_leadjet++;
            ptSumTransverse2_leadjet += pT;
          }

          if (fabs(dphi_leadjet) < PI / 3.) {  // Toward region
            nTowards_leadjet++;
            ptSumTowards_leadjet += pT;
          }

          if (fabs(dphi_leadjet) > 2. * PI / 3.) {  // Away region
            nAway_leadjet++;
            ptSumAway_leadjet += pT;
          }

        }  //jet found

        if (!p_leadtrack.isZero()) {
          double dphi_leadtrack = p.phi() - phileadtrack;
          while (dphi_leadtrack > PI) {
            dphi_leadtrack = dphi_leadtrack - 2.0 * PI;
          }
          while (dphi_leadtrack < -PI) {
            dphi_leadtrack = dphi_leadtrack + 2. * PI;
          }

          if (dphi_leadtrack > PI / 3. && dphi_leadtrack < PI * 2. / 3.) {  // Transverse1 region
            nTransverse_leadtrack++;
            ptSumTransverse_leadtrack += pT;
            nTransverse1_leadtrack++;
            ptSumTransverse1_leadtrack += pT;
          }

          if (dphi_leadtrack < -PI / 3. && dphi_leadtrack > -PI * 2. / 3.) {  // Transverse2 region
            nTransverse_leadtrack++;
            ptSumTransverse_leadtrack += pT;
            nTransverse2_leadtrack++;
            ptSumTransverse2_leadtrack += pT;
          }

          if (fabs(dphi_leadtrack) < PI / 3.) {  // Toward region
            nTowards_leadtrack++;
            ptSumTowards_leadtrack += pT;
          }

          if (fabs(dphi_leadtrack) > 2. * PI / 3.) {  // Away region
            nAway_leadtrack++;
            ptSumAway_leadtrack += pT;
          }

        }  //track found

      }  //Loop over particles

      const double fullarea = 8. / 3. * PI;
      const double halfarea = 4. / 3. * PI;

      if (!p_leadjet.isZero()) {
        if (nTransverse2_leadjet > nTransverse1_leadjet) {
          nTransverseMax_leadjet = nTransverse2_leadjet;
          nTransverseMin_leadjet = nTransverse1_leadjet;
        }

        else {
          nTransverseMax_leadjet = nTransverse1_leadjet;
          nTransverseMin_leadjet = nTransverse2_leadjet;
        }

        if (ptSumTransverse2_leadjet > ptSumTransverse1_leadjet) {
          ptSumTransverseMax_leadjet = ptSumTransverse2_leadjet;
          ptSumTransverseMin_leadjet = ptSumTransverse1_leadjet;
        }

        else {
          ptSumTransverseMax_leadjet = ptSumTransverse1_leadjet;
          ptSumTransverseMin_leadjet = ptSumTransverse2_leadjet;
        }

        _h_Nch_vs_leadJetPt_transDiff->fill(pTleadjet / GeV,
                                            1. / halfarea * (nTransverseMax_leadjet - nTransverseMin_leadjet));
        _h_PtSum_vs_leadJetPt_transDiff->fill(
            pTleadjet / GeV, 1. / halfarea * (ptSumTransverseMax_leadjet - ptSumTransverseMin_leadjet));
        _h_Nch_vs_leadJetPt_transAvg->fill(pTleadjet / GeV,
                                           1. / fullarea * (nTransverseMax_leadjet + nTransverseMin_leadjet));
        _h_PtSum_vs_leadJetPt_transAvg->fill(pTleadjet / GeV,
                                             1. / fullarea * (ptSumTransverseMax_leadjet + ptSumTransverseMin_leadjet));
        _h_Nch_vs_leadJetPt_transMax->fill(pTleadjet / GeV, 1. / halfarea * nTransverseMax_leadjet);
        _h_PtSum_vs_leadJetPt_transMax->fill(pTleadjet / GeV, 1. / halfarea * ptSumTransverseMax_leadjet);
        _h_Nch_vs_leadJetPt_transMin->fill(pTleadjet / GeV, 1. / halfarea * nTransverseMin_leadjet);
        _h_PtSum_vs_leadJetPt_transMin->fill(pTleadjet / GeV, 1. / halfarea * ptSumTransverseMin_leadjet);

      }  //for leading jet

      if (!p_leadtrack.isZero()) {
        if (nTransverse2_leadtrack > nTransverse1_leadtrack) {
          nTransverseMax_leadtrack = nTransverse2_leadtrack;
          nTransverseMin_leadtrack = nTransverse1_leadtrack;
        }

        else {
          nTransverseMax_leadtrack = nTransverse1_leadtrack;
          nTransverseMin_leadtrack = nTransverse2_leadtrack;
        }

        if (ptSumTransverse2_leadtrack > ptSumTransverse1_leadtrack) {
          ptSumTransverseMax_leadtrack = ptSumTransverse2_leadtrack;
          ptSumTransverseMin_leadtrack = ptSumTransverse1_leadtrack;
        }

        else {
          ptSumTransverseMax_leadtrack = ptSumTransverse1_leadtrack;
          ptSumTransverseMin_leadtrack = ptSumTransverse2_leadtrack;
        }

        _h_Nch_vs_leadTrackPt_transDiff->fill(pTleadtrack / GeV,
                                              1. / halfarea * (nTransverseMax_leadtrack - nTransverseMin_leadtrack));
        _h_PtSum_vs_leadTrackPt_transDiff->fill(
            pTleadtrack / GeV, 1. / halfarea * (ptSumTransverseMax_leadtrack - ptSumTransverseMin_leadtrack));
        _h_Nch_vs_leadTrackPt_transAvg->fill(pTleadtrack / GeV,
                                             1. / fullarea * (nTransverseMax_leadtrack + nTransverseMin_leadtrack));
        _h_PtSum_vs_leadTrackPt_transAvg->fill(
            pTleadtrack / GeV, 1. / fullarea * (ptSumTransverseMax_leadtrack + ptSumTransverseMin_leadtrack));

        _h_Nch_vs_leadTrackPt_transMax->fill(pTleadtrack / GeV, 1. / halfarea * nTransverseMax_leadtrack);
        _h_PtSum_vs_leadTrackPt_transMax->fill(pTleadtrack / GeV, 1. / halfarea * ptSumTransverseMax_leadtrack);
        _h_Nch_vs_leadTrackPt_transMin->fill(pTleadtrack / GeV, 1. / halfarea * nTransverseMin_leadtrack);
        _h_PtSum_vs_leadTrackPt_transMin->fill(pTleadtrack / GeV, 1. / halfarea * ptSumTransverseMin_leadtrack);

      }  //for leading track
    }

    /// Normalise histograms etc., after the run
    void finalize() {}

  private:
    Profile1DPtr _h_Nch_vs_leadJetPt_transMax;
    Profile1DPtr _h_PtSum_vs_leadJetPt_transMax;
    Profile1DPtr _h_Nch_vs_leadJetPt_transMin;
    Profile1DPtr _h_PtSum_vs_leadJetPt_transMin;
    Profile1DPtr _h_Nch_vs_leadJetPt_transDiff;
    Profile1DPtr _h_PtSum_vs_leadJetPt_transDiff;
    Profile1DPtr _h_Nch_vs_leadJetPt_transAvg;
    Profile1DPtr _h_PtSum_vs_leadJetPt_transAvg;

    Profile1DPtr _h_Nch_vs_leadTrackPt_transMax;
    Profile1DPtr _h_PtSum_vs_leadTrackPt_transMax;
    Profile1DPtr _h_Nch_vs_leadTrackPt_transMin;
    Profile1DPtr _h_PtSum_vs_leadTrackPt_transMin;
    Profile1DPtr _h_Nch_vs_leadTrackPt_transDiff;
    Profile1DPtr _h_PtSum_vs_leadTrackPt_transDiff;
    Profile1DPtr _h_Nch_vs_leadTrackPt_transAvg;
    Profile1DPtr _h_PtSum_vs_leadTrackPt_transAvg;
  };

  // This global object acts as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2015_PAS_FSQ_15_007);
}  // namespace Rivet
