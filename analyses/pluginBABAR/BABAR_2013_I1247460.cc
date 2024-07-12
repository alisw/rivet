// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B+ -> omega l+ nu_l
  class BABAR_2013_I1247460 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2013_I1247460);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==PID::BPLUS), "UFS");

      // Book histograms
      book(_h_q2 ,1, 1, 1);
      book(_nB,"TMP/nB");
    }
    
    // Calculate the Q2 using mother and daughter charged lepton
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
      // Get B+ Mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	_nB->fill();
        if (isSemileptonicDecay(p, {PID::OMEGA, PID::POSITRON, PID::NU_E}) ||
            isSemileptonicDecay(p, {PID::OMEGA, PID::ANTIMUON, PID::NU_MU})) {
	  _h_q2->fill(q2(p,PID::OMEGA));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // BR in units of 10^{-5}
      scale(_h_q2,1e5/ *_nB);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nB;
    Histo1DPtr _h_q2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2013_I1247460);

}
