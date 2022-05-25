// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Search for new physics with dijet angular distributions at 13 TeV
  class CMS_2017_I1519995 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1519995);


    /// Book projections and histograms
    void init() {
      FastJets antikt(FinalState(), FastJets::ANTIKT, 0.4);
      declare(antikt, "ANTIKT");
      /// @todo We need a better way!!
      {Histo1DPtr tmp; _h_chi_dijet.add(4800., 8000., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(4200., 4800., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3600., 4200., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3000., 3600., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(2400., 3000., book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1900., 2400., book(tmp, 6, 1, 1));}
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const FourMomentum j0(jets[0].mom()), j1(jets[1].mom());
      if (fabs(j0.rap()+j1.rap())/2 > 1.11) vetoEvent;

      const double mjj = (j0+j1).mass();
      const double chi = exp(fabs(j0.rap()-j1.rap()));
      if (chi < 16) _h_chi_dijet.fill(mjj/GeV, chi, 1.0);
    }


    /// Normalize histograms
    void finalize() {
      for (Histo1DPtr hist : _h_chi_dijet.histos()) normalize(hist);
    }


    BinnedHistogram _h_chi_dijet;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2017_I1519995);

}
