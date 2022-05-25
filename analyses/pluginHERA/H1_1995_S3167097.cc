// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/CentralEtHCM.hh"

namespace Rivet {


  /// H1 energy flow in DIS
  ///
  /// @todo Make histograms match those in HepData and use autobooking
  ///
  /// @author Leif Lonnblad
  /// @author Andy Buckley
  class H1_1995_S3167097 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1995_S3167097);


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      const DISKinematics& diskin = declare(DISKinematics(), "Kinematics");
      const DISFinalState& fshcm = declare(DISFinalState(diskin, DISFinalState::BoostFrame::HCM), "FS");
      declare(CentralEtHCM(fshcm), "Y1HCM");

      // Histograms
      /// @todo Convert to use autobooking and correspond to HepData data tables

      _hEtFlow.resize(9);
      for (size_t i = 0; i < 9; ++i) {
        book(_sumw[i], "sumW_" + to_str(i));
        book(_hEtFlow[i] ,to_str(i), 24, -6, 6);
      }
      book(_tmphAvEt, "TMP/hAvEt", 9, 1.0, 10.0);
      book(_tmphAvX , "TMP/hAvX", 9, 1.0, 10.0);
      book(_tmphAvQ2, "TMP/hAvQ2", 9, 1.0, 10.0);
      book(_tmphN   , "TMP/hN", 9, 1.0, 10.0);
    }


    /// Calculate the bin number from the DISKinematics projection
    /// @todo Convert to use a HEPUtils Binning1D
    size_t _getbin(const DISKinematics& dk) {
      if (inRange(dk.Q2()/GeV2, 5.0, 10.0)) {
        if (inRange(dk.x(), 1e-4, 2e-4)) return 0;
        if (inRange(dk.x(), 2e-4, 5e-4) && dk.Q2() > 6.0*GeV2) return 1;
      } else if (inRange(dk.Q2()/GeV2, 10.0, 20.0)) {
        if (inRange(dk.x(), 2e-4, 5e-4)) return 2;
        if (inRange(dk.x(), 5e-4, 8e-4)) return 3;
        if (inRange(dk.x(), 8e-4, 1.5e-3)) return 4;
        if (inRange(dk.x(), 1.5e-3, 4e-3)) return 5;
      } else if (inRange(dk.Q2()/GeV2, 20.0, 50.0)) {
        if (inRange(dk.x(), 5e-4, 1.4e-3)) return 6;
        if (inRange(dk.x(), 1.4e-3, 3e-3)) return 7;
        if (inRange(dk.x(), 3e-3, 1e-2)) return 8;
      }
      return -1;
    }


    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) vetoEvent;
      const CentralEtHCM& y1 = apply<CentralEtHCM>(event, "Y1HCM");
      if ( y1.failed() ) vetoEvent;

      const int ibin = _getbin(dk);
      if (ibin < 0) vetoEvent;

      _sumw[ibin]->fill();

      for (size_t i = 0, N = fs.particles().size(); i < N; ++i) {
        const double rap = fs.particles()[i].rapidity();
        const double et = fs.particles()[i].Et();
        _hEtFlow[ibin]->fill(rap, et/GeV);
      }

      /// @todo Use fillBin?
      _tmphAvEt->fill(ibin + 1.5, y1.sumEt()/GeV);
      _tmphAvX->fill(ibin + 1.5, dk.x());
      _tmphAvQ2->fill(ibin + 1.5, dk.Q2()/GeV2);
      _tmphN->fill(ibin + 1.5);
    }


    void finalize() {
      for (size_t ibin = 0; ibin < 9; ++ibin)
        scale(_hEtFlow[ibin], 0.5/ *_sumw[ibin]);
      /// @todo Improve this!
      Scatter2DPtr s21,s22,s23;
      divide(_tmphAvEt,_tmphN,s21);
      book(s21, "21");
      divide(_tmphAvX,_tmphN,s22);
      book(s22, "22");
      divide(_tmphAvQ2,_tmphN,s23);
      book(s23, "23");
    }

    //@}


  private:

    /// Histograms for the \f$ E_T \f$ flow
    vector<Histo1DPtr> _hEtFlow;

    /// Temporary histograms for averages in different kinematical bins.
    Histo1DPtr _tmphAvEt, _tmphAvX, _tmphAvQ2, _tmphN;

    /// Weights counters for each kinematic bin
    array<CounterPtr, 9> _sumw;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(H1_1995_S3167097, H1_1995_I396365);

}
