// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TASSO_1990_S2148048 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1990_S2148048);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::pT >=  0.1/GeV);
      declare(cfs, "CFS");

      // Thrust and sphericity
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Histos
      int offset = 0;
      switch (int(sqrtS()/GeV)) {
        case 14:
          offset = 0;
          break;
        case 22:
          offset = 1;
          break;
        case 35:
          offset = 2;
          break;
        case 44:
          offset = 3;
          break;
      }
      book(_h_xp[0]      , 2, 1, 1+offset);
      book(_h_xp[1]      , 3, 1, 1+offset);
      book(_h_xi         , 4, 1, 1+offset);
      book(_h_pT         , 5, 1, 1+offset);
      book(_h_sphericity , 6, 1, 1+offset);
      book(_h_aplanarity , 7, 1, 1+offset);
      book(_h_thrust     , 8, 1, 1+offset);
      book(_sumWPassed,"/TMP/_sumWPassed");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // TASSO hadronic event selection TODO: move this into a trigger definition
      // See page 2 in publication
      // Condition 1)  --- require at least 5 (4) 'good' tracks
      int nch = cfs.particles().size();
      if ( (int(sqrtS()/GeV) > 27 && nch < 5) || (int(sqrtS()/GeV) <= 27 && nch < 4 ) ) {
        MSG_DEBUG("Failed # good tracks cut: " << nch);
        vetoEvent;
      }
      // Condition 2) ---
      // Condition 5) --- scalar momentum (not pT!!!) sum >= 0.265*s
      double momsum = 0.0;
      for (const Particle& p : cfs.particles()) {
        const double mom = p.p3().mod();
        momsum += mom;
      }
      if (momsum <=0.265 * sqrtS()/GeV) {
        MSG_DEBUG("Failed pTsum cut: " << momsum << " < " << 0.265 * sqrtS()/GeV);
        vetoEvent;
      }

      // Raise counter for events that pass trigger conditions
      _sumWPassed->fill();

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      //const Vector3 & thrustAxis = thrust.thrustAxis ();
      //double theta = thrustAxis.theta();
      //if ( fabs(cos(theta)) >= 0.8 ) {
        //MSG_DEBUG("Failed thrust angle cut: " << fabs(cos(theta)));
        //vetoEvent;
      //}

      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      // Fill histograms in order of appearance in paper
      for (const Particle& p : cfs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        // Scaled momenta.
        const double mom = mom3.mod();
        const double scaledMom = 2.*mom/sqrtS();
        const double pTin = dot(mom3, sphericity.sphericityMajorAxis());
        const double pTout = dot(mom3, sphericity.sphericityMinorAxis());
	const double pT=sqrt(sqr(pTin)+sqr(pTout));
        _h_xp[0]->fill(scaledMom);
        _h_xp[1]->fill(scaledMom);
	_h_xi   ->fill(-log(scaledMom));
	_h_pT   ->fill(pT);
      }
      // event shapes
      _h_sphericity->fill(sphericity.sphericity());
      _h_aplanarity->fill(sphericity.aplanarity());
      _h_thrust->fill(thrust.thrust());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_xp[0], 1./ *_sumWPassed);
      scale(_h_xp[1], 1./ *_sumWPassed);
      scale(_h_xi   , 1./ *_sumWPassed);
      scale(_h_pT   , 1./ *_sumWPassed);
      normalize(_h_sphericity);
      normalize(_h_aplanarity);
      normalize(_h_thrust    );
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_xp[2], _h_xi, _h_pT, _h_sphericity, _h_aplanarity, _h_thrust;
    CounterPtr _sumWPassed;
    //@}


  };



  RIVET_DECLARE_ALIASED_PLUGIN(TASSO_1990_S2148048, TASSO_1990_I294755);

}
