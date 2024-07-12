// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 dijet invariant mass measurement
  class D0_2010_S8566488 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_S8566488);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      FastJets conefinder(fs, FastJets::D0ILCONE, 0.7);
      declare(conefinder, "ConeFinder");

      {Histo1DPtr tmp; _h_m_dijet.add(0.0, 0.4, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_m_dijet.add(0.4, 0.8, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_m_dijet.add(0.8, 1.2, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_m_dijet.add(1.2, 1.6, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_m_dijet.add(1.6, 2.0, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_m_dijet.add(2.0, 2.4, book(tmp, 6, 1, 1));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      const Jets& jets = apply<JetAlg>(e, "ConeFinder").jetsByPt(40.0*GeV);
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double ymax = std::max(j0.absrap(), j1.absrap());
      double mjj = FourMomentum(j0+j1).mass();

      _h_m_dijet.fill(ymax, mjj/TeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      _h_m_dijet.scale(crossSection()/sumOfWeights(), this);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    BinnedHistogram _h_m_dijet;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_S8566488, D0_2010_I846483);

}
