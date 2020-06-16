// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief q^2 in D0 -> K- e+ nu_e decays
  class BABAR_2007_I1091435 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2007_I1091435);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_q2 ,1, 1, 1);
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
      // Loop over B0 mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::D0)) {
        if (isSemileptonicDecay(p, {PID::KMINUS, PID::POSITRON, PID::NU_E}) ) {
	  _h_q2->fill(q2(p, PID::KMINUS));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_q2, 1.);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_q2;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2007_I1091435);

}
