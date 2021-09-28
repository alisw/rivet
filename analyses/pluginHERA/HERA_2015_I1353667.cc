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
      book(_h_pTD, 1, 1, 1);
      book(_h_etaD, 2, 1, 1);
      book(_h_zD, 3, 1, 1);
      book(_h_Q2, 4, 1, 1);
      book(_h_y, 5, 1, 1);
 book(_h_Q2y, 6, 1, 1);
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
        _h_pTD->fill(dstar.pT()/GeV);
        _h_etaD->fill(dstar.eta());
        _h_zD->fill(zD/GeV);
        _h_Q2->fill(kin.Q2()/GeV2);
        _h_y->fill(kin.y());
      }

      // // Double-differential (y,Q2) histograms
      // _h_Q2y->fill(kin.Q2()/GeV2, kin.y());

    }


    /// Normalise histograms etc., after the run
    void finalize() { 
      const double sf = crossSection()/nanobarn/sumOfWeights();
      scale(_h_pTD, sf);
      scale(_h_etaD, sf);
      scale(_h_zD, sf);
      scale(_h_Q2, sf);
      scale(_h_y, sf);
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
