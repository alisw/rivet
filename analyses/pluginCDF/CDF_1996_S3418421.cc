// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF dijet angular distributions
  class CDF_1996_S3418421 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1996_S3418421);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");

      {Histo1DPtr tmp; _h_chi.add(241.0, 300.0, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_chi.add(300.0, 400.0, book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_chi.add(400.0, 517.0, book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_chi.add(517.0, 625.0, book(tmp, 1, 1, 4));}
      {Histo1DPtr tmp; _h_chi.add(625.0,1800.0, book(tmp, 1, 1, 5));}
      book(_h_ratio,  2, 1, 1);
      book(_htmp_chi_above_25 ,"TMP/chiabove25", refData(2, 1, 1));
      book(_htmp_chi_below_25 ,"TMP/chibelow25", refData(2, 1, 1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(50.0*GeV);
      if (jets.size() < 2) vetoEvent;

      const FourMomentum jet1 = jets[0].momentum();
      const FourMomentum jet2 = jets[1].momentum();
      const double eta1 = jet1.eta();
      const double eta2 = jet2.eta();
      const double chi = exp(fabs(eta1 - eta2));
      if (fabs(eta2) > 2.0 || fabs(eta1) > 2.0 || chi > 5.0) vetoEvent;

      double m = FourMomentum(jet1 + jet2).mass();
      _h_chi.fill(m, chi, weight);

      // Fill ratio numerator or denominator depending on chi value
      ((chi > 2.5) ? _htmp_chi_above_25 : _htmp_chi_below_25)->fill(m/GeV, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr hist : _h_chi.histos()) {
        normalize(hist);
      }
      divide(_htmp_chi_below_25, _htmp_chi_above_25, _h_ratio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    BinnedHistogram _h_chi;
    Histo1DPtr _htmp_chi_above_25, _htmp_chi_below_25;
    Scatter2DPtr _h_ratio;
    //@}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_1996_S3418421, CDF_1996_I423414);

}
