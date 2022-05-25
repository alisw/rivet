// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISFinalState.hh"

namespace Rivet {


/// @brief Measurement of $D^{\pm}$ and $D^0$ production in deep inelastic scattering using a lifetime tag at HERA
/// @author Andrii Verbytskyi
class ZEUS_2008_I810112 : public Analysis {

  public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2008_I810112);

      /// @name Analysis methods
      ///@{

      /// Book histograms and initialise projections before the run
      void init() {
          /// Initialise and register projections
          FinalState fs;
          const DISKinematics& diskin = DISKinematics();
          declare(diskin,"Kinematics");
          declare(UnstableParticles(), "UPS");


          book(_h_Dp_q2, 3, 1, 1);
          book(_h_Dp_x, 4, 1, 1);
          book(_h_Dp_pt, 5, 1, 1);
          book(_h_Dp_eta, 6, 1, 1);

          book(_h_Dp_yinq2[0], "TMP/Dp0", refData( 11, 1, 1));
          book(_h_Dp_yinq2[1], "TMP/Dp1", refData( 11, 1, 2));
          book(_h_Dp_yinq2[2], "TMP/Dp2", refData( 11, 1, 3));

          book(_s_Dp_yinq2[0], 11, 1, 1);
          book(_s_Dp_yinq2[1], 11, 1, 2);
          book(_s_Dp_yinq2[2], 11, 1, 3);

          book(_h_D0_q2, 7, 1, 1);
          book(_h_D0_x, 8, 1, 1);
          book(_h_D0_pt, 9, 1, 1);
          book(_h_D0_eta, 10, 1, 1);

          book(_h_D0_yinq2[0], "TMP/D00", refData( 12, 1, 1));
          book(_h_D0_yinq2[1], "TMP/D01", refData( 12, 1, 2));
          book(_h_D0_yinq2[2], "TMP/D02", refData( 12, 1, 3));

          book(_s_D0_yinq2[0], 12, 1, 1);
          book(_s_D0_yinq2[1], 12, 1, 2);
          book(_s_D0_yinq2[2], 12, 1, 3);


      }
  /// A routine that selects the bin in kinematic space
      int _getbinQ2_OK(const DISKinematics& dk) {
          if (inRange(dk.Q2()/GeV2, 5.0, 9.0)) return 0;
          if (inRange(dk.Q2()/GeV2, 9.0, 44.0)) return 1;
          if (inRange(dk.Q2()/GeV2, 44.0, 1000.0)) return 2;
          return -1;
      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {
          /// DIS kinematics
          const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
          double q2  = dk.Q2();
          double x   = dk.x();
          double y   = dk.y();

          ///Assure the event falls into the kinematic range of the measurement.
          if (!inRange(q2/GeV2, 5, 1000)) vetoEvent;
          if (!inRange(y, 0.02, 0.7)) vetoEvent;

          /// Find out in which bin the measurement falls.
          int bin = _getbinQ2_OK(dk);
          if ( bin < 0 ) vetoEvent;

          const UnstableParticles& ufs = apply<UnstableParticles>(event, "UPS");

          /// Get \f$D^0$\f particles
          for (const Particle& p : filter_select(ufs.particles(), Cuts::abspid == PID::D0)) {
              ///But not not from \f$D^{*\pm}$\f decays
              if (p.hasAncestor(PID::DSTARPLUS)) continue;
              if (p.hasAncestor(PID::DSTARMINUS)) continue;
              /// Select particles only in the \f$\eta-p_{T}$\f region
              if (!inRange(p.eta(), -1.6, 1.6)) continue;
              if (!inRange(p.pt()/GeV, 1.5, 15.0)) continue;
              /// Fill the general histograms
              _h_D0_pt->fill(p.pt()/GeV);
              _h_D0_eta->fill(p.eta());
              _h_D0_q2->fill(q2/GeV2);
              _h_D0_x->fill(x);
              /// Fill the cross-sections in selected kinematic bins
              _h_D0_yinq2[bin]->fill(y);
          }

          /// Get \f$D^{\pm}$\f particles
          for (const Particle& p : filter_select(ufs.particles(), Cuts::abspid == PID::DPLUS)) {
              /// Select particles only in the \f$\eta-p_{T}$\f region
              if (!inRange(p.eta(), -1.6, 1.6)) continue;
              if (!inRange(p.pt()/GeV, 1.5, 15.0)) continue;
              /// Fill the general histograms
              _h_Dp_pt->fill(p.pt()/GeV);
              _h_Dp_eta->fill(p.eta());
              _h_Dp_q2->fill(q2/GeV2);
              _h_Dp_x->fill(x);
              /// Fill the cross-sections in selected kinematic bins
              _h_Dp_yinq2[bin]->fill(y);
          }
      }
      /// Normalise histograms etc., after the run
      void finalize() {
          const double sf = crossSection()/nanobarn/sumOfWeights();
          /// For CS comparison the MC should be corrected to difference between BR used for data and in MC
          scale(_h_Dp_pt, sf);
          scale(_h_Dp_eta, sf);
          scale(_h_Dp_q2, sf);
          scale(_h_Dp_x, sf);

          for (size_t i = 0; i < 3; i++) {
              scale(_h_Dp_yinq2[i], sf);
              barchart(_h_Dp_yinq2[i], _s_Dp_yinq2[i], false);
          }

          scale(_h_D0_pt, sf);
          scale(_h_D0_eta, sf);
          scale(_h_D0_q2, sf);
          scale(_h_D0_x, sf);
          
          for (size_t i = 0; i < 3; i++) {
              scale(_h_D0_yinq2[i], sf);
              barchart(_h_D0_yinq2[i], _s_D0_yinq2[i], false);
          }
      }
      ///@}

    private:

      Histo1DPtr
      _h_D0_pt,
      _h_D0_eta,
      _h_D0_q2,
      _h_D0_x,
      _h_D0_yinq2[3];
      Histo1DPtr
      _h_Dp_pt,
      _h_Dp_eta,
      _h_Dp_q2,
      _h_Dp_x,
      _h_Dp_yinq2[3];
      
      Scatter2DPtr _s_Dp_yinq2[3];
      Scatter2DPtr _s_D0_yinq2[3];

  };

  /// The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ZEUS_2008_I810112);
}
