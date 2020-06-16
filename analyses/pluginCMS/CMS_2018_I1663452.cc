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
      declare(antikt, "ANTIKT");
      {Histo1DPtr tmp; _h_chi_dijet.add(6000., 13000., book(tmp,1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(5400., 6000., book(tmp,2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(4800., 5400., book(tmp,3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(4200., 4800., book(tmp,4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3600., 4200., book(tmp,5, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(3000., 3600., book(tmp,6, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(2400., 3000., book(tmp,7, 1, 1));}
    }

    void analyze(const Event& event) {
      const Jets& jets = applyProjection<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;
      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0+y1)/2. > 1.11) vetoEvent;
      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet.fill(mjj, chi);
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
  DECLARE_RIVET_PLUGIN(CMS_2018_I1663452);

}
