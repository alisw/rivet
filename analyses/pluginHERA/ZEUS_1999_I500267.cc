// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

/// @brief  Measurement of High-$Q^2$ Neutral-Current in DIS
class ZEUS_1999_I500267 : public Analysis {

  public:
      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1999_I500267);

      /// @name Analysis methods
      ///@{

      /// Book histograms and initialise projections before the run
      void init() {

          const DISKinematics& diskin = DISKinematics();
          declare(diskin,"Kinematics");

          book(_h_Q2,   1, 1, 1);
          book(_h_x[0], 2, 1, 1);
          book(_h_x[1], 3, 1, 1);
          book(_h_x[2], 4, 1, 1);

      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {

          /// DIS kinematics
          const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");

          double q2  = dk.Q2();
          double x   = dk.x();
          double y   = dk.y();

          if (y > 0.95) vetoEvent;
          if (q2 < 400) vetoEvent;
          if (q2 > 400) _h_x[0]->fill(x);
          if (q2 > 2500) _h_x[1]->fill(x);
          if (q2 > 10000) _h_x[2]->fill(x);
          _h_Q2->fill(q2);

      }

      /// Normalise histograms etc., after the run
      void finalize() {
          const double norm = crossSection()/picobarn/sumOfWeights();
          scale(_h_Q2, norm);
          scale(_h_x[0],norm);
          scale(_h_x[1], norm);
          scale(_h_x[2], norm);
      }

    private:

      Histo1DPtr _h_Q2;
      Histo1DPtr _h_x[3];
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ZEUS_1999_I500267);


}
