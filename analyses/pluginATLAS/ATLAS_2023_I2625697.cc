// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief TEEC in multijets at 13 TeV
  class ATLAS_2023_I2625697 : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2625697);

    /// @name Analysis methods
    //@{

    void init() {
      // Projections AntiKt4TruthJets
      const FinalState fs(Cuts::abseta < 4.5);
      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "Jets");

      // Book histograms
      for (int k = 0; k < 11; ++k) {
        book(_hEEC[k], k+1, 1, 1);
        book(_hAEEC[k], k+12, 1, 1);
      }
    }

    void analyze(const Event& event) {

      const double htBins[11] = {1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2300.0, 2600.0, 3000.0, 3500.0, 13000.0};

      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 60.0*GeV && Cuts::abseta < 2.4);

      if (jets.size() < 2) vetoEvent;

      const double sumPt12 = jets[0].pt()+jets[1].pt();
      if (sumPt12 < 1000.0*GeV) vetoEvent;

      const double sumEt = sum(jets, Kin::Et, 0.0);

      for (const Jet& j1 : jets) {
        for (const Jet& j2 : jets) {

          const double etWeight = j1.Et()*j2.Et()/sqr(sumEt);
          const double dPhi = deltaPhi(j1, j2);

          double cosPhi = cos(dPhi);
          if (cos(dPhi) == 1.0) cosPhi = 0.9999;

          _hEEC[0]->fill(cosPhi, etWeight);

          for (int k = 0; k < 10; ++k){
            if (inRange(sumPt12/GeV, htBins[k], htBins[k+1])) _hEEC[k+1]->fill(cosPhi, etWeight);
          }
        }
      }
    }

    void finalize(){

      normalize(_hEEC);

      for (int iBin = 0; iBin < 11; ++iBin) {
        size_t nPoints = _hEEC[iBin]->numBins();
        for (size_t k = 0; k < nPoints/2; ++k) {
          const double x = _hEEC[iBin]->bin(k).midpoint();
          const double ex = _hEEC[iBin]->bin(k).xWidth()/2;
          const double y = _hEEC[iBin]->bin(k).height() - _hEEC[iBin]->bin(nPoints-(k+1)).height();
          const double ey = sqrt( sqr(_hEEC[iBin]->bin(k).heightErr()) + sqr(_hEEC[iBin]->bin(nPoints-(k+1)).heightErr()) );
          _hAEEC[iBin]->addPoint(x, y, ex, ey);
        }
       }
    }

    //@}

  private:

    map<int, Histo1DPtr> _hEEC;
    map<int, Scatter2DPtr> _hAEEC;
  };

  // Declare the plugin
  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2625697);

}
