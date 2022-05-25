// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Measurement of charged particle transverse momentum spectra in deep inelastic scattering
  class H1_1997_I424463 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1997_I424463);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      DISLepton dl;
      declare(dl, "Lepton");
      DISKinematics diskin;
      declare(diskin, "Kinematics");
      ChargedFinalState cfs(dl.remainingFinalState());
      declare(cfs, "CFS");

      for (int offset = 0; offset < 10; offset++) {
        book(_h_pt[offset], 1 + offset, 1, 1);
        book(_h_eta[offset], 29 + offset, 1, 1);
        book(_h_pt[offset+10], 11 + offset, 1, 1);
      }
      book(_h_Q2, "TMP/Q2", 10, 0, 10);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const double q2  = dk.Q2();
      const double x   = dk.x();
      const double y   = dk.y();

      if (!inRange(q2/GeV2, 5.0, 50.0)) vetoEvent;
      if (!inRange(y, 0.05, 1.0)) vetoEvent;
      if (!inRange(x, 0.0001, 0.01)) vetoEvent;

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.particles().size() < 2)  vetoEvent;

      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      const Vector4 leptonMom = dl.out().momentum();
      const double enel = leptonMom.t();
      const double thel = leptonMom.angle(dk.beamHadron().mom())/degree;

      bool cut = enel/GeV > 12. && thel > 157. && thel < 173;
      if (!cut) vetoEvent;

      const int of = _getbinQ2x(dk);
      if ( of < 0 ) vetoEvent;
      const Vector3 gMom = (dl.out().momentum() - dl.in().momentum()).p3().unit();
      if ( of > 0 ) _h_Q2->fill(of);
      _h_Q2->fill(0);

      for (const auto  p: cfs.particles()) {
        const double thp = p.momentum().angle(dk.beamHadron().mom())/degree;

        if (!( thp > 8. && thp < 155)) continue;

        const double costhetastar = dot(p.momentum().p3().unit(),gMom);
        const double etastar = -0.5*log(( 1.0 - costhetastar)/(1.0 + costhetastar));
        if (inRange(etastar, 0.5, 2.5)) {
          int delta = -1;
          if (inRange(etastar, 0.5, 1.5)) delta = 0;
          if (inRange(etastar, 1.5, 2.5)) delta = 1;
          if (delta > -1) {
            if (of > 0) {
              _h_pt[of + 10*delta]->fill(p.pt()/GeV);
            }
            _h_pt[10*delta]->fill(p.pt()/GeV);
          }
        }
        if (p.pt()/GeV > 1.0) {
          const double eta = p.momentum().p3().eta();//-0.5*log((1-costheta)/(1+costheta));
          if ( of > 0 ) _h_eta[of]->fill(eta);
          _h_eta[0]->fill(eta);
        }

      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int offset = 0; offset < 10; offset++) {
        const double st = _h_Q2->binAt(offset).height();
        scale(_h_eta[offset], 1.0/std::max(1.0, st));
        scale(_h_pt[offset], 1.0/std::max(1.0, st));
        scale(_h_pt[offset+10], 1.0/std::max(1.0, st));
      }
    }

    /// @}


    /// Utility function to get appropriate (Q2,x) bin indices
    int _getbinQ2x(const DISKinematics& dk) {
      const double q2 = dk.Q2()/GeV2, x = dk.x();
      if (inRange(q2, 5.0, 10.0) && inRange(x, 0.0001, 0.0002)) return 1;
      if (inRange(q2, 6.0, 10.0) && inRange(x, 0.0002, 0.0005)) return 2;
      if (inRange(q2, 10.0, 20.0)) {
        if (inRange(x, 0.0002, 0.0005)) return 3;
        if (inRange(x, 0.0005, 0.0008)) return 4;
        if (inRange(x, 0.0008, 0.0015)) return 5;
        if (inRange(x, 0.0015, 0.0040)) return 6;
      }
      if (inRange(q2, 20.0, 50.0)) {
        if (inRange(x, 0.0005, 0.0014)) return 7;
        if (inRange(x, 0.0014, 0.0030)) return 8;
        if (inRange(x, 0.0030, 0.010)) return 9;
      }
      if (inRange(q2, 5.0, 50.0) && inRange(x, 0.0001, 0.01) ) return 0;
      return -1;
    }


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pt[20];
    Histo1DPtr _h_eta[10];
    Histo1DPtr _h_Q2;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(H1_1997_I424463);

}
