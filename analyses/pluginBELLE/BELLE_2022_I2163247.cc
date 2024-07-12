// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  B0 -> pi - l+ nu_l
  class BELLE_2022_I2163247 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2163247);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==PID::B0), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) book(_h[ix],1,1,1+ix);
      book(_nB,"TMP/nb");
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
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	_nB->fill();
        if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E})) {
	  double qq = q2(p, PID::PIMINUS);
	  _h[0]->fill(qq);
	  _h[2]->fill(qq);
	}
	else if(isSemileptonicDecay(p, {PID::PIMINUS, PID::ANTIMUON, PID::NU_MU})) {
	  double qq = q2(p, PID::PIMINUS);
	  _h[1]->fill(qq);
	  _h[2]->fill(qq);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	double fact = 1e4;
	if(ix==2) fact /=2;
	scale(_h[ix], fact/ *_nB);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    CounterPtr _nB;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2163247);

}
