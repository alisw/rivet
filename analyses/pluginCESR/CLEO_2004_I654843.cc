// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D -> pi,K semi-leptonic q^2
  class CLEO_2004_I654843 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2004_I654843);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==PID::D0), "UFS");
      // histograms
      book(_h_q2_D0_K ,1,1,1);
      book(_h_q2_D0_pi,1,1,2);
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::abspid==abs(mesonID))[0];
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
      // Loop over D mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	if(isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E}) ||
	   isSemileptonicDecay(p, {PID::PIPLUS , PID::ELECTRON, PID::NU_EBAR}) )
	  _h_q2_D0_pi->fill(q2(p, PID::PIMINUS));
	else if(isSemileptonicDecay(p, {PID::KMINUS, PID::POSITRON, PID::NU_E}) ||
		isSemileptonicDecay(p, {PID::KPLUS , PID::ELECTRON, PID::NU_EBAR}))
	  _h_q2_D0_K ->fill(q2(p, PID::KMINUS));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_q2_D0_K );
      normalize(_h_q2_D0_pi);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr  _h_q2_D0_K, _h_q2_D0_pi;
    CounterPtr _nD0,_nDp;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2004_I654843);

}
