// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Dijet angular distributions and search for quark compositeness at 7 TeV
  class CMS_2011_S8968497 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8968497);


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      {Histo1DPtr tmp; _h_chi_dijet.add(2200., 7000., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1800., 2200., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1400., 1800., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1100., 1400., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 850., 1100., book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 650.,  850., book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 500.,  650., book(tmp, 7, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 350.,  500., book(tmp, 8, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add( 250.,  350., book(tmp, 9, 1, 1));}
    }


    void analyze(const Event& event) {
      const double weight = 1.0;
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
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
      for (Histo1DPtr hist : _h_chi_dijet.histos()) {
        normalize(hist);
      }
    }


  private:

    BinnedHistogram _h_chi_dijet;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8968497, CMS_2011_I889175);

}
