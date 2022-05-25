// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS UE charged particles vs. leading jet at 2.76 TeV
  class CMS_2015_I1385107 : public Analysis {
  public:
    /// Constructor
    CMS_2015_I1385107() : Analysis("CMS_2015_I1385107"),
                          ETACUT(2.0),
                          AREATOT(2*ETACUT * 2*M_PI),
                          AREA3(AREATOT / 3.),
                          AREA6(AREATOT / 6.)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(Cuts::abseta < 2 && Cuts::pT > 500*MeV);
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV);
      const FastJets jetpro(cfsforjet, FastJets::SISCONE, 0.5);
      declare(jetpro, "Jets");

      book(_h_Nch_TransAVE_vs_pT ,1, 1, 1); // Nch vs. pT_max      (TransAVE)
      book(_h_Sum_TransAVE_vs_pT ,2, 1, 1); // sum(pT) vs. pT_max  (TransAVE)
      book(_h_Nch_TransMAX_vs_pT ,3, 1, 1); // Nch vs. pT_max      (TransMAX)
      book(_h_Sum_TransMAX_vs_pT ,4, 1, 1); // sum(pT) vs. pT_max  (TransMAX)
      book(_h_Nch_TransMIN_vs_pT ,5, 1, 1); // Nch vs. pT_max      (TransMIN)
      book(_h_Sum_TransMIN_vs_pT ,6, 1, 1); // sum(pT) vs. pT_max  (TransMIN)
      book(_h_Nch_TransDIF_vs_pT ,7, 1, 1); // Nch vs. pT_max      (TransDIF)
      book(_h_Sum_TransDIF_vs_pT ,8, 1, 1); // sum(pT) vs. pT_max  (TransDIF)
    }


    /// Local definition of a signed dphi, for use in differentating L and R trans regions
    double signedDeltaPhi(double jetphi, double partphi) {
      double delta = partphi - jetphi;
      while (delta <= -PI) delta += 2 * PI;
      while (delta > PI) delta -= 2 * PI;
      return delta;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_lead;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(1*GeV)) {
        if (j.abseta() < 2.0) {
          p_lead = j.momentum();
          break;
        }
      }
      if (p_lead.isZero()) vetoEvent;
      const double phi_lead = p_lead.phi();
      const double pT_lead  = p_lead.pT();

      // Loop on charged particles and separate Left and Right transverse regions
      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();
      int nch_TransLeft = 0, nch_TransRight = 0;
      double ptSum_TransLeft = 0., ptSum_TransRight = 0.;
      for (const Particle& p : particles) {
        const double dphi = signedDeltaPhi(phi_lead, p.momentum().phi());
        if (!inRange(fabs(dphi), PI/3, 2*PI/3.)) continue; //< only fill trans regions
        if (dphi < 0) {  // Transverse Right region
          nch_TransRight += 1;
          ptSum_TransRight += p.pT() / GeV;
        } else if (dphi > 0) {  // Transverse Left region
          nch_TransLeft += 1;
          ptSum_TransLeft += p.pT() / GeV;
        }
      }

      // Translate to min and max (+sum and diff) Transverse regions
      const int nch_TransMIN = std::min(nch_TransLeft, nch_TransRight);
      const int nch_TransMAX = std::max(nch_TransLeft, nch_TransRight);
      const int nch_TransSUM = nch_TransMAX + nch_TransMIN;
      const int nch_TransDIF = nch_TransMAX - nch_TransMIN;
      //
      const double ptSum_TransMIN = std::min(ptSum_TransLeft, ptSum_TransRight);
      const double ptSum_TransMAX = std::max(ptSum_TransLeft, ptSum_TransRight);
      const double ptSum_TransSUM = ptSum_TransMAX + ptSum_TransMIN;
      const double ptSum_TransDIF = ptSum_TransMAX - ptSum_TransMIN;

      // Fill profiles
      _h_Nch_TransMIN_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransMIN);
      _h_Sum_TransMIN_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransMIN);
      //
      _h_Nch_TransMAX_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransMAX);
      _h_Sum_TransMAX_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransMAX);
      //
      _h_Nch_TransAVE_vs_pT->fill(pT_lead/GeV, 1/AREA3 * nch_TransSUM);
      _h_Sum_TransAVE_vs_pT->fill(pT_lead/GeV, 1/AREA3 * ptSum_TransSUM);
      //
      _h_Nch_TransDIF_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransDIF);
      _h_Sum_TransDIF_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransDIF);
    }


  private:

    // Data members like post-cuts event weight counters go here
    const double ETACUT, AREATOT, AREA3, AREA6;

    /// Histograms
    Profile1DPtr _h_Nch_TransAVE_vs_pT, _h_Sum_TransAVE_vs_pT;
    Profile1DPtr _h_Nch_TransDIF_vs_pT, _h_Sum_TransDIF_vs_pT;
    Profile1DPtr _h_Nch_TransMIN_vs_pT, _h_Sum_TransMIN_vs_pT;
    Profile1DPtr _h_Nch_TransMAX_vs_pT, _h_Sum_TransMAX_vs_pT;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2015_I1385107);

}
