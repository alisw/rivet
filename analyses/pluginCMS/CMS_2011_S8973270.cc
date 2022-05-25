// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// B-Bbar angular correlations based on secondary vertex reconstruction
  class CMS_2011_S8973270 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8973270);


    /// @name Analysis methods
    /// @{

    void init() {
      FinalState fs;
      FastJets jetproj(fs, FastJets::ANTIKT, 0.5);
      jetproj.useInvisibles();
      declare(jetproj, "Jets");

      UnstableParticles ufs;
      declare(ufs, "UFS");

      // Book histograms
      book(_h_dsigma_dR_56GeV ,1,1,1);
      book(_h_dsigma_dR_84GeV ,2,1,1);
      book(_h_dsigma_dR_120GeV ,3,1,1);
      book(_h_dsigma_dPhi_56GeV ,4,1,1);
      book(_h_dsigma_dPhi_84GeV ,5,1,1);
      book(_h_dsigma_dPhi_120GeV ,6,1,1);

      book(_c["MCDR56"],     "_MCDR56");
      book(_c["MCDR84"],     "_MCDR84");
      book(_c["MCDR120"],    "_MCDR120");
      book(_c["MCDPhi56"],   "_MCDPhi56");
      book(_c["MCDPhi84"],   "_MCDPhi84");
      book(_c["MCDPhi120"], "_MCDPhi120");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      const Jets& jets = apply<FastJets>(event,"Jets").jetsByPt();
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      // Find the leading jet pT and eta
      if (jets.size() == 0) vetoEvent;
      const double ljpT = jets[0].pT();
      const double ljeta = jets[0].eta();
      MSG_DEBUG("Leading jet pT / eta: " << ljpT << " / " << ljeta);

      // Minimum requirement for event
      if (ljpT > 56*GeV && fabs(ljeta) < 3.0) {
        // Find B hadrons in event
        int nab = 0, nb = 0; //counters for all B and independent B hadrons
        double etaB1 = 7.7, etaB2 = 7.7;
        double phiB1 = 7.7, phiB2 = 7.7;
        double pTB1 = 7.7, pTB2 = 7.7;

        for (const Particle& p : ufs.particles()) {
          int aid = p.abspid();
          if (aid/100 == 5 || aid/1000==5) {
            nab++;
            // 2J+1 == 1 (mesons) or 2 (baryons)
            if (aid%10 == 1 || aid%10 == 2) {
              // No B decaying to B
              if (aid != 5222 && aid != 5112 && aid != 5212 && aid != 5322) {
                if (nb==0) {
                  etaB1 = p.eta();
                  phiB1 = p.phi();
                  pTB1 = p.pT();
                } else if (nb==1) {
                  etaB2 = p.eta();
                  phiB2 = p.phi();
                  pTB2 = p.pT();
                }
                nb++;
              }
            }
            MSG_DEBUG("ID " << aid <<  " B hadron");
          }
        }

        if (nb==2 && pTB1 > 15*GeV && pTB2 > 15*GeV && fabs(etaB1) < 2.0 && fabs(etaB2) < 2.0) {
          double dPhi = deltaPhi(phiB1, phiB2);
          double dR = deltaR(etaB1, phiB1, etaB2, phiB2);
          MSG_DEBUG("DR/DPhi " << dR << " " << dPhi);

          // MC counters
          if (dR > 2.4) _c["MCDR56"]->fill();
          if (dR > 2.4 && ljpT > 84*GeV) _c["MCDR84"]->fill();
          if (dR > 2.4 && ljpT > 120*GeV) _c["MCDR120"]->fill();
          if (dPhi > 3.*PI/4.) _c["MCDPhi56"]->fill();
          if (dPhi > 3.*PI/4. && ljpT > 84*GeV) _c["MCDPhi84"]->fill();
          if (dPhi > 3.*PI/4. && ljpT > 120*GeV) _c["MCDPhi120"]->fill();

          _h_dsigma_dR_56GeV->fill(dR, weight);
          if (ljpT > 84*GeV) _h_dsigma_dR_84GeV->fill(dR, weight);
          if (ljpT > 120*GeV) _h_dsigma_dR_120GeV->fill(dR, weight);
          _h_dsigma_dPhi_56GeV->fill(dPhi, weight);
          if (ljpT > 84*GeV) _h_dsigma_dPhi_84GeV->fill(dPhi, weight);
          if (ljpT > 120*GeV) _h_dsigma_dPhi_120GeV->fill(dPhi, weight);
          //MSG_DEBUG("nb " << nb << " " << nab);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG("crossSection " << crossSection() << " sumOfWeights " << sumOfWeights());

      // Hardcoded bin widths
      double DRbin = 0.4;
      double DPhibin = PI/8.0;
      // Find out the correct numbers
      double nDataDR56 = 25862.20;
      double nDataDR84 = 5675.55;
      double nDataDR120 = 1042.72;
      double nDataDPhi56 = 24220.00;
      double nDataDPhi84 = 4964.00;
      double nDataDPhi120 = 919.10;
      double normDR56 = safediv(nDataDR56, dbl(*_c["MCDR56"]), crossSection()/sumOfWeights());
      double normDR84 = safediv(nDataDR84, dbl(*_c["MCDR84"]), crossSection()/sumOfWeights());
      double normDR120 = safediv(nDataDR120, dbl(*_c["MCDR120"]), crossSection()/sumOfWeights());
      double normDPhi56 = safediv(nDataDPhi56, dbl(*_c["MCDPhi56"]), crossSection()/sumOfWeights());
      double normDPhi84 = safediv(nDataDPhi84, dbl(*_c["MCDPhi84"]), crossSection()/sumOfWeights());
      double normDPhi120 = safediv(nDataDPhi120, dbl(*_c["MCDPhi120"]), crossSection()/sumOfWeights());
      scale(_h_dsigma_dR_56GeV, normDR56*DRbin);
      scale(_h_dsigma_dR_84GeV, normDR84*DRbin);
      scale(_h_dsigma_dR_120GeV, normDR120*DRbin);
      scale(_h_dsigma_dPhi_56GeV, normDPhi56*DPhibin);
      scale(_h_dsigma_dPhi_84GeV, normDPhi84*DPhibin);
      scale(_h_dsigma_dPhi_120GeV, normDPhi120*DPhibin);
    }

    /// @}


  private:

    /// Counters
    map<string, CounterPtr> _c;

    /// @name Histograms
    /// @{
    Histo1DPtr _h_dsigma_dR_56GeV, _h_dsigma_dR_84GeV, _h_dsigma_dR_120GeV;
    Histo1DPtr _h_dsigma_dPhi_56GeV, _h_dsigma_dPhi_84GeV, _h_dsigma_dPhi_120GeV;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8973270, CMS_2011_I889807);

}
