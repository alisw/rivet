// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2015_I1327224 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1327224);


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      {Histo1DPtr tmp; _h_chi_dijet.add(4200., 8000., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3600., 4200., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3000., 3600., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(2400., 3000., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1900., 2400., book(tmp, 5, 1, 1));}
    }


    void analyze(const Event& event) {
      const double weight = 1.0;
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0 + y1) / 2. > 1.11) vetoEvent;

      double mjj = FourMomentum(j0 + j1).mass();
      if (mjj/GeV <1900) vetoEvent;

      double chi = exp(fabs(y0 - y1));
      if (chi >= 16.) vetoEvent;

      // Fill the histogram
      _h_chi_dijet.fill(mjj/GeV, chi, weight);
    }

    void finalize() {
      for (Histo1DPtr hist : _h_chi_dijet.histos()) {
        normalize(hist);
      }
    }

  private:
    BinnedHistogram _h_chi_dijet;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2015_I1327224);
}
