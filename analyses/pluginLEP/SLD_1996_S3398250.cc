// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include <cmath>

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief SLD multiplicities at mZ
  /// @author Peter Richardson
  class SLD_1996_S3398250 : public Analysis {
  public:

    /// Constructor
    SLD_1996_S3398250()
      : Analysis("SLD_1996_S3398250")
    {}

    /// @name Analysis methods
    //@{


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");

      book(_h_bottom ,1, 1, 1);
      book(_h_charm  ,2, 1, 1);
      book(_h_light  ,3, 1, 1);

      book(_weightLight, "_weightLight");
      book(_weightCharm, "_weightCharm");
      book(_weightBottom, "_weightBottom");

      book(scatter_c, 4,1,1);
      book(scatter_b, 5,1,1);

    }


    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1: case 2: case 3:
        _weightLight ->fill();
        _h_light->fillBin(0, numParticles);
        break;
      case 4:
        _weightCharm ->fill();
        _h_charm->fillBin(0, numParticles);
        break;
      case 5:
        _weightBottom->fill();
        _h_bottom->fillBin(0, numParticles);
        break;
      }

    }


    void multiplicity_subtract(const Histo1DPtr first, const Histo1DPtr second, Scatter2DPtr & scatter) {
      const double x  = first->bin(0).xMid();
      const double ex = first->bin(0).xWidth()/2.;
      const double y  = first->bin(0).area() - second->bin(0).area();
      const double ey = sqrt(sqr(first->bin(0).areaErr()) + sqr(second->bin(0).areaErr()));
      scatter->addPoint(x, y, ex, ey);
    }


    void finalize() {
      if (_weightBottom->val() != 0) scale(_h_bottom, 1./ *_weightBottom);
      if (_weightCharm->val()  != 0) scale(_h_charm,  1./ *_weightCharm );
      if (_weightLight->val()  != 0) scale(_h_light,  1./ *_weightLight );

      multiplicity_subtract(_h_charm,  _h_light, scatter_c);
      multiplicity_subtract(_h_bottom, _h_light, scatter_b);
    }

    //@}


  private:

    Scatter2DPtr scatter_c, scatter_b;
    /// @name Weights
    //@{
    CounterPtr _weightLight;
    CounterPtr _weightCharm;
    CounterPtr _weightBottom;
    //@}

    Histo1DPtr _h_bottom;
    Histo1DPtr _h_charm;
    Histo1DPtr _h_light;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_1996_S3398250);

}
