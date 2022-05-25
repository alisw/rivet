// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 dijet angular distributions
  class D0_2009_S8320160 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2009_S8320160);


    /// @name Analysis methods
    /// @{

    // Book histograms
    void init() {
      FinalState fs;
      FastJets conefinder(fs, FastJets::D0ILCONE, 0.7);
      declare(conefinder, "ConeFinder");

      {Histo1DPtr tmp; _h_chi_dijet.add(250., 300., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(300., 400., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(400., 500., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(500., 600., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(600., 700., book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(700., 800., book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(800., 900., book(tmp, 7, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(900.,1000., book(tmp, 8, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1000.,1100.,book(tmp, 9, 1, 1));}
      {Histo1DPtr tmp; _h_chi_dijet.add(1100.,1960, book(tmp, 10, 1, 1));}
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = 1.0;

      const Jets& jets = apply<JetAlg>(e, "ConeFinder").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();

      if (fabs(y0+y1)>2) vetoEvent;

      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet.fill(mjj, chi, weight);
    }


    /// Finalize
    void finalize() {
      for (Histo1DPtr hist : _h_chi_dijet.histos()) {
        normalize(hist);
      }
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    BinnedHistogram _h_chi_dijet;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2009_S8320160, D0_2009_I824127);

}
