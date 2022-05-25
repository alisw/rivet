// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Based on BABAR_2015_I1334693
  class BESIII_2015_I1391138 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1391138);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2_K, 1, 1, 3);
      book(_h_q2_pi, 2, 1, 3);
      book(nD0, "TMP/DCounter");
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }
    
    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over D0 mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::D0)) {
	nD0->fill();
        if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E})) {
          _h_q2_pi->fill(q2(p, PID::PIMINUS));
        }
	else if(isSemileptonicDecay(p, {PID::KMINUS, PID::POSITRON, PID::NU_E})) {
          _h_q2_K ->fill(q2(p, PID::KMINUS));
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // scale by D0 lifetime = 410.1e-6 ps (from PDG 2014 used in paper)
      // and bin width 0.1 K and 0.2 pi
      scale(_h_q2_K , 1./nD0->sumW()/410.1e-6*0.1);
      scale(_h_q2_pi, 1./nD0->sumW()/410.1e-6*0.2);
    }
    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_q2_K,_h_q2_pi;
    CounterPtr nD0;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2015_I1391138);


}
