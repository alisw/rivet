// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2012_I1090423 : public Analysis {
  public:

    CMS_2012_I1090423()
      : Analysis("CMS_2012_I1090423")
    { }


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      {Histo1DPtr tmp; _h_chi_dijet.add(3000, 7000, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(2400, 3000, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1900, 2400, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1500, 1900, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1200, 1500, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1000, 1200, book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 800, 1000, book(tmp, 7, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 600,  800, book(tmp, 8, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 400,  600, book(tmp, 9, 1, 1));}
    }


    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const double y0 = jets[0].rapidity();
      const double y1 = jets[1].rapidity();
      if (fabs(y0+y1)/2 > 1.11) vetoEvent;

      const double chi = exp(fabs(y0-y1));
      if (chi > 16) vetoEvent;

      const FourMomentum jj = jets[0].momentum() + jets[1].momentum();
       _h_chi_dijet.fill(jj.mass(), chi, 1.0);
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
  DECLARE_RIVET_PLUGIN(CMS_2012_I1090423);

}
