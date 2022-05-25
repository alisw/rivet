// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Inelastic xsection in pp collisions at 13 TeV for charged particles in LHCb acceptance
  class LHCB_2018_I1665223 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2018_I1665223);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Register projection
      declare(ChargedFinalState(Cuts::etaIn(ETAMIN, ETAMAX)), "lbCFS");

      // Book histogram
      book(_h_ppInel, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Eliminate non-inelastic events and empty events in LHCb
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "lbCFS");
      if (cfs.particles().size() == 0) vetoEvent;

      for (const Particle& myp : cfs.particles()) {
        if (hasLongLivedParent(myp)) continue;
        if (!isLongLivedParticle(myp)) continue;
        // Only particles with p > 2 GeV measured
        if (myp.momentum().p() < PMIN) continue;
        // Histo gets filled only for inelastic events (at least one prompt charged particle)
        _h_ppInel->fill(sqrtS());
        break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ppInel, crossSection()/millibarn/sumOfWeights()); // norm to cross section (one-sided LHCb)
    }

    //@}


  	bool isLongLivedParticle(const Particle& p) {
  		// Stable long-lived final-state charged particles in LHCb
  		static const int stablePids[9] = {11,13,211,321,2212,3112,3222,3312,3334};
  		for (int stablePid : stablePids) {
          if (p.abspid() == stablePid) return true;
  		}
  		return false;
  	}

    bool hasLongLivedParent(const Particle& p) {
      // List of PDG IDs for particle with lifetimes higher than 0.03 ns (3.E-11 s) - long lived particles according to LHC MB&UE WG
      static const int longLivedPids[20] = {3334,3322,3312,3222,3122,3112,2212,2112,321,310,130,211,20022,480000000,11,12,13,14,16,22};
      for (int longLivedPid : longLivedPids) {
        if (p.hasParentWith(Cuts::abspid == longLivedPid)) return true;
      }
      return false;
    }

    /// @name Histogram
    Histo1DPtr _h_ppInel;

    /// Cut constants
    const double ETAMIN = 2.0, ETAMAX = 5.0, PMIN = 2.0*GeV;

  };


  RIVET_DECLARE_PLUGIN(LHCB_2018_I1665223);

}
