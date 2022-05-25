// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross-sections of \f$K_{0}$\f and \f$\Lambda$\f in DIS
  class H1_2009_I810046 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(H1_2009_I810046);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const DISKinematics& diskin = DISKinematics();
      declare(diskin, "Kinematics");
      declare(UnstableParticles(), "UPS");

      book(_h_K0S_q2, 4, 1, 1);
      book(_h_K0S_x, 5, 1, 1);
      book(_h_K0S_pt, 6, 1, 1);
      book(_h_K0S_eta, 7, 1, 1);

      book(_h_LAMBDA_q2, 8, 1, 1);
      book(_h_LAMBDA_x, 9, 1, 1);
      book(_h_LAMBDA_pt, 10, 1, 1);
      book(_h_LAMBDA_eta, 11, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const double q2  = dk.Q2();
      const double x   = dk.x();
      const double y   = dk.y();
      const int orientation = dk.orientation();

      if (!inRange(q2/GeV2, 2.0, 100.0)) vetoEvent;
      if (!inRange(y, 0.1, 0.6)) vetoEvent;
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UPS");

      for (const Particle& p: filter_select(ufs.particles(), Cuts::abspid == abs(PID::K0S))) {
        if (!inRange(p.pt()/GeV, 0.5, 3.5)) continue;
        if (!inRange(p.eta(), -1.3, 1.3)) continue;
        _h_K0S_q2->fill(q2/GeV2);
        _h_K0S_x->fill(x);
        _h_K0S_pt->fill(p.pt()/GeV);
        _h_K0S_eta->fill(p.eta()*orientation);
      }

      for (const Particle& p: filter_select(ufs.particles(), Cuts::abspid == abs(PID::LAMBDA))) {
        if (!inRange(p.pt()/GeV, 0.5, 3.5)) continue;
        if (!inRange(p.eta(), -1.3, 1.3)) continue;
        _h_LAMBDA_q2->fill(q2/GeV2);
        _h_LAMBDA_x->fill(x);
        _h_LAMBDA_pt->fill(p.pt()/GeV);
        _h_LAMBDA_eta->fill(p.eta()*orientation);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/nanobarn/sumOfWeights();
      scale( _h_K0S_pt, sf);
      scale( _h_K0S_eta, sf);
      scale( _h_K0S_q2, sf);
      scale( _h_K0S_x, sf/1000);

      scale( _h_LAMBDA_pt, sf);
      scale( _h_LAMBDA_eta, sf);
      scale( _h_LAMBDA_q2, sf);
      scale( _h_LAMBDA_x, sf/1000);
    }

    /// @}

    /// @name Histograms
    /// @}
    Histo1DPtr _h_K0S_pt, _h_K0S_eta, _h_K0S_x, _h_K0S_q2;
    Histo1DPtr _h_LAMBDA_pt, _h_LAMBDA_eta, _h_LAMBDA_x, _h_LAMBDA_q2;
    /// @}

  };


  DECLARE_RIVET_PLUGIN(H1_2009_I810046);

}
