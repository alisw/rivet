// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// CMS azimuthal decorrelations
  class CMS_2011_S8950903 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8950903);


    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.5);
      declare(akt, "antikT");

      {Histo1DPtr tmp; _h_deltaPhi.add( 80.,  110., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.add(110.,  140., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.add(140.,  200., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.add(200.,  300., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.add(300., 7000., book(tmp, 5, 1, 1));}
    }


    void analyze(const Event & event) {
      const double weight = 1.0;

      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      if (fabs(jets[0].eta()) > 1.1 || jets[0].pT() < 80.) vetoEvent;
      if (fabs(jets[1].eta()) > 1.1 || jets[1].pT() < 30.) vetoEvent;

      double dphi = deltaPhi(jets[0].momentum(), jets[1].phi());

      _h_deltaPhi.fill(jets[0].pT(), dphi, weight);
    }


    void finalize() {
      normalize(_h_deltaPhi.histos(), 1.);
    }


  private:

    BinnedHistogram _h_deltaPhi;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8950903, CMS_2011_I885663);

}
