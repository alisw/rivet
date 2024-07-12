// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// UE charged particles vs. leading jet
  class CMS_2011_S9120041 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S9120041);


    void init() {
      const ChargedFinalState cfs((Cuts::etaIn(-2.0, 2.0) && Cuts::pT >=  500*MeV));
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  500*MeV));
      const FastJets jetpro(cfsforjet, FastJets::SISCONE, 0.5);
      declare(jetpro, "Jets");

      if (isCompatibleWithSqrtS(7000.)) {
        book(_h_Nch_vs_pT ,1, 1, 1); // Nch vs. pT_max
        book(_h_Sum_vs_pT ,2, 1, 1); // sum(pT) vs. pT_max
        book(_h_pT3_Nch   ,5, 1, 1);   // transverse Nch,     pT_max > 3GeV
        book(_h_pT3_Sum   ,6, 1, 1);   // transverse sum(pT), pT_max > 3GeV
        book(_h_pT3_pT    ,7, 1, 1);   // transverse pT,      pT_max > 3GeV
        book(_h_pT20_Nch  ,8, 1, 1);   // transverse Nch,     pT_max > 20GeV
        book(_h_pT20_Sum  ,9, 1, 1);   // transverse sum(pT), pT_max > 20GeV
        book(_h_pT20_pT   ,10, 1, 1);  // transverse pT,      pT_max > 20GeV
      }

      if (isCompatibleWithSqrtS(900.)) {
        book(_h_Nch_vs_pT ,3, 1, 1); // Nch vs. pT_max
        book(_h_Sum_vs_pT ,4, 1, 1); // sum(pT) vs. pT_max
        book(_h_pT3_Nch   ,11, 1, 1);  // transverse Nch,     pT_max > 3GeV
        book(_h_pT3_Sum   ,12, 1, 1);  // transverse sum(pT), pT_max > 3GeV
        book(_h_pT3_pT    ,13, 1, 1);  // transverse pT,      pT_max > 3GeV
      }

      book(sumOfWeights3, "TMP/sumOfWeights3");
      book(sumOfWeights20, "TMP/sumOfWeights20");
      book(_nch_tot_pT3, "TMP/nch_tot_pT3");
      book(_nch_tot_pT20, "TMP/nch_tot_pT20");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_lead;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(1.0*GeV)) {
        if (j.abseta() < 2.0) {
          p_lead = j.momentum();
          break;
        }
      }
      if (p_lead.isZero()) vetoEvent;
      const double philead = p_lead.phi();
      const double pTlead  = p_lead.pT();

      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();

      int nTransverse = 0;
      double ptSumTransverse = 0.;
      for (const Particle& p : particles) {
        double dphi = fabs(deltaPhi(philead, p.phi()));
        if (dphi>PI/3. && dphi<PI*2./3.) {   // Transverse region
          nTransverse++;

          const double pT = p.pT()/GeV;
          ptSumTransverse += pT;

          if (pTlead > 3.0*GeV) _h_pT3_pT->fill(pT);
          if (isCompatibleWithSqrtS(7000.) && pTlead > 20.0*GeV) _h_pT20_pT->fill(pT);
        }
      }

      const double area = 8./3. * PI;
      _h_Nch_vs_pT->fill(pTlead/GeV, 1./area*nTransverse);
      _h_Sum_vs_pT->fill(pTlead/GeV, 1./area*ptSumTransverse);
      if(pTlead > 3.0*GeV) {
        _h_pT3_Nch->fill(nTransverse);
        _h_pT3_Sum->fill(ptSumTransverse);
        sumOfWeights3->fill();
        _nch_tot_pT3->fill(nTransverse);
      }
      if (isCompatibleWithSqrtS(7000.) && pTlead > 20.0*GeV) {
        _h_pT20_Nch->fill(nTransverse);
        _h_pT20_Sum->fill(ptSumTransverse);
        sumOfWeights20->fill();
        _nch_tot_pT20->fill(nTransverse);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pT3_Nch);
      normalize(_h_pT3_Sum);
      if (sumOfWeights3->val() != 0.0) normalize(_h_pT3_pT, *_nch_tot_pT3 / *sumOfWeights3);

      if (isCompatibleWithSqrtS(7000.)) {
        normalize(_h_pT20_Nch);
        normalize(_h_pT20_Sum);
        if (sumOfWeights20->val() != 0.0) normalize(_h_pT20_pT, *_nch_tot_pT20 / *sumOfWeights20);
      }
    }



  private:

    /// @{
    CounterPtr sumOfWeights3;
    CounterPtr sumOfWeights20;

    CounterPtr _nch_tot_pT3;
    CounterPtr _nch_tot_pT20;

    Profile1DPtr _h_Nch_vs_pT;
    Profile1DPtr _h_Sum_vs_pT;
    Histo1DPtr _h_pT3_Nch;
    Histo1DPtr _h_pT3_Sum;
    Histo1DPtr _h_pT3_pT;
    Histo1DPtr _h_pT20_Nch;
    Histo1DPtr _h_pT20_Sum;
    Histo1DPtr _h_pT20_pT;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S9120041, CMS_2011_I916908);

}
