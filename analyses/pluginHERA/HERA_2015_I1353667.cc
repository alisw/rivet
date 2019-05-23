// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Combined H1/ZEUS D* production cross-sections in DIS
  class HERA_2015_I1353667 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(HERA_2015_I1353667);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      // FinalState fs;
      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "Dstars");
      //Cuts::abspid == PID::DSTARPLUS

      // Book histograms
      _h_pTD = bookHisto1D(1, 1, 1);
      _h_etaD = bookHisto1D(2, 1, 1);
      _h_zD = bookHisto1D(3, 1, 1);
      _h_Q2 = bookHisto1D(4, 1, 1);
      _h_y = bookHisto1D(5, 1, 1);
      // _h_Q2y = bookHisto2D(6, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Determine kinematics, including event orientation
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");
      //const int orientation = kin.orientation();

      // Q2 and inelasticity cuts
      if (!inRange(kin.Q2(), 1.5*GeV2, 1000*GeV2)) vetoEvent;
      if (!inRange(kin.y(), 0.02, 0.7)) vetoEvent;


      // D* reconstruction
      const Particles unstables = apply<ParticleFinder>(event, "Dstars")
        .particles(Cuts::pT > 1.5*GeV && Cuts::abseta < 1.5);
      const Particles dstars = filter_select(unstables, [](const Particle& p){ return p.abspid() == PID::DSTARPLUS; });
      if (dstars.empty()) vetoEvent;
      MSG_DEBUG("#D* = " << dstars.size());
      const Particle& dstar = dstars.front();
      const double zD = (dstar.E() - dstar.pz()) / (2*kin.beamLepton().E()*kin.y());

      // Single-differential histograms with higher low-Q2 cut
      if (kin.Q2() > 5*GeV2) {
        _h_pTD->fill(dstar.pT()/GeV, event.weight());
        _h_etaD->fill(dstar.eta(), event.weight());
        _h_zD->fill(zD/GeV, event.weight());
        _h_Q2->fill(kin.Q2()/GeV2, event.weight());
        _h_y->fill(kin.y(), event.weight());
      }

      // // Double-differential (y,Q2) histograms
      // _h_Q2y->fill(kin.Q2()/GeV2, kin.y(), event.weight());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale({_h_pTD, _h_etaD, _h_zD, _h_Q2, _h_y}, crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pTD, _h_etaD, _h_zD, _h_Q2, _h_y;
    Histo2DPtr _h_Q2y;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HERA_2015_I1353667);


}
