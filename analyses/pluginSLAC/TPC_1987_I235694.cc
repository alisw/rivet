// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief TPC flavour separated N charged at 29 GeV
  class TPC_1987_I235694 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TPC_1987_I235694);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Histograms
      book(_h_all   , 5, 1, 4);
      book(_h_light , 4, 1, 4);
      book(_h_charm , 3, 1, 4);
      book(_h_bottom, 2, 1, 4);
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");
    }


    /// Perform the per-event analysis
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
	_h_light->fill(sqrtS()/GeV,numParticles);
        break;
      case 4:
	_h_charm->fill(sqrtS()/GeV,numParticles);
        break;
      case 5:
	_h_bottom->fill(sqrtS()/GeV,numParticles);
        break;
      }
      _h_all->fill(sqrtS()/GeV,numParticles);
    }


    /// Normalise histograms etc., after the run
    void finalize() { }
    //@}


  private:

    /// @name Multiplicities
    //@{
    Profile1DPtr _h_all;
    Profile1DPtr _h_light;
    Profile1DPtr _h_charm;
    Profile1DPtr _h_bottom;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TPC_1987_I235694);


}
