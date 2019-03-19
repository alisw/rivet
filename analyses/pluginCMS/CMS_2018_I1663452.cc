// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2018_I1663452 : public Analysis {
  public:

    CMS_2018_I1663452()
      : Analysis("CMS_2018_I1663452")
    { }


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.4);
      addProjection(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(6000., 13000., bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(5400., 6000., bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(4800., 5400., bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(4200., 4800., bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram(3600., 4200., bookHisto1D(5, 1, 1));
      _h_chi_dijet.addHistogram(3000., 3600., bookHisto1D(6, 1, 1));
      _h_chi_dijet.addHistogram(2400., 3000., bookHisto1D(7, 1, 1));
    }

    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets& jets = applyProjection<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;
      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0+y1)/2. > 1.11) vetoEvent;
      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet.fill(mjj, chi, weight);
    }


    void finalize() {
      foreach (Histo1DPtr hist, _h_chi_dijet.getHistograms()) {
        normalize(hist);
      }
    }


  private:

    BinnedHistogram<double> _h_chi_dijet;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_I1663452);

}
