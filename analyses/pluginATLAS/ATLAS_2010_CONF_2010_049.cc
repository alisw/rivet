// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2010_CONF_2010_049 : public Analysis {
  public:

    ATLAS_2010_CONF_2010_049()
      : Analysis("ATLAS_2010_CONF_2010_049")
    {    }


    void init() {
      ChargedFinalState cfs((Cuts::etaIn(-1.5, 1.5) && Cuts::pT >=  0.5*GeV));
      declare(cfs, "CFS");

      FastJets jetsproj6(cfs, FastJets::ANTIKT, 0.6);
      declare(jetsproj6, "Jets6");

      FastJets jetsproj4(cfs, FastJets::ANTIKT, 0.4);
      declare(jetsproj4, "Jets4");

      // @todo tmp YOs
      for (size_t i=0 ; i<2 ; i++) {
        book(_h_xsec[i]       ,1+i, 1, 1);
        book(_h_frag_04_06[i] ,3+i, 1, 1);
        book(_h_frag_06_10[i] ,3+i, 2, 1);
        book(_h_frag_10_15[i] ,3+i, 3, 1);
        book(_h_frag_15_24[i] ,3+i, 4, 1);
        book(_njets_04_06[i], "njets_04_06_"+to_string(i));
        book(_njets_06_10[i], "njets_06_10_"+to_string(i));
        book(_njets_10_15[i], "njets_10_15_"+to_string(i));
        book(_njets_15_24[i], "njets_15_24_"+to_string(i));
      }
    }


    void analyze(const Event& event) {
      const FastJets & jetsproj6 = apply<FastJets>(event, "Jets6");
      const FastJets & jetsproj4 = apply<FastJets>(event, "Jets4");
      Jets alljets[2];
      alljets[0] = jetsproj6.jetsByPt(4.0*GeV);
      alljets[1] = jetsproj4.jetsByPt(4.0*GeV);

      for (size_t i=0 ; i<2 ; i++) {
        Jets jets;

        // First we want to make sure that we only use jets within |eta|<0.57
        for (const Jet& jet : alljets[i]) {
          if (jet.abseta()<0.57) {
            jets.push_back(jet);
          }
        }
        for (const Jet& jet : jets) {
          const double pTjet = jet.pT();
          const double pjet = jet.p3().mod();
          _h_xsec[i]->fill(pTjet);
          if (pTjet > 24*GeV) continue;
          for (const Particle& p : jet.particles()) {
            double z = p.p3().mod()/pjet;
            if (z >= 1) z = 0.9999; // Make sure that z=1 doesn't go into overflow
            if (pTjet > 15*GeV) {
              _h_frag_15_24[i]->fill(z);
            }
            else if (pTjet > 10*GeV) {
              _h_frag_10_15[i]->fill(z);
            }
            else if (pTjet > 6*GeV) {
              _h_frag_06_10[i]->fill(z);
            }
            else {
              _h_frag_04_06[i]->fill(z);
            }
          }
          if (pTjet > 15*GeV) {
            _njets_15_24[i]->fill();
          }
          else if (pTjet > 10*GeV) {
            _njets_10_15[i]->fill();
          }
          else if (pTjet > 6*GeV) {
            _njets_06_10[i]->fill();
          }
          else {
            _njets_04_06[i]->fill();
          }
        }
      }
    }

    void finalize() {
      for (size_t i=0 ; i<2 ; i++) {
        // deta = 2*0.57
        scale(_h_xsec[i], crossSection()/microbarn/sumOfWeights()/(2*0.57));
        scale(_h_frag_04_06[i], 1./_njets_04_06[i]->val());
        scale(_h_frag_06_10[i], 1./_njets_06_10[i]->val());
        scale(_h_frag_10_15[i], 1./_njets_10_15[i]->val());
        scale(_h_frag_15_24[i], 1./_njets_15_24[i]->val());
      }
    }


  private:

    Histo1DPtr _h_xsec[2];
    Histo1DPtr _h_frag_04_06[2];
    Histo1DPtr _h_frag_06_10[2];
    Histo1DPtr _h_frag_10_15[2];
    Histo1DPtr _h_frag_15_24[2];
    CounterPtr _njets_04_06[2];
    CounterPtr _njets_06_10[2];
    CounterPtr _njets_10_15[2];
    CounterPtr _njets_15_24[2];
  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2010_CONF_2010_049);

}
