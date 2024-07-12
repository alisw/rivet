// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> D semileptonic
  class BELLE_2022_I2169621 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2169621);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(Cuts::pid==511 or
      				Cuts::pid==521), "UFS");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
      	book(_h[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<2;++ix)
      	book(_nB[ix],"TMP/nB_"+toString(ix));
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }
    
    // Calculate the recoil w using mother and daugher meson
    double recoilW(const Particle& B, int mesonID) {
      // TODO why does that not work with const?
      Particle D = filter_select(B.children(), Cuts::pid==mesonID)[0];
      FourMomentum q = B.mom() - D.mom();
      return (B.mom()*B.mom() + D.mom()*D.mom() - q*q )/ (2. * sqrt(B.mom()*B.mom()) * sqrt(D.mom()*D.mom()) );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles()) {
      	if(p.children().size()<=1) continue;
      	if(p.pid()==PID::BPLUS) {
      	  _nB[0]->fill();
      	  if (isSemileptonicDecay(p, {PID::D0BAR,PID::POSITRON,PID::NU_E}))  _h[0]->fill(recoilW(p, PID::D0BAR));
      	  if (isSemileptonicDecay(p, {PID::D0BAR,PID::ANTIMUON,PID::NU_MU})) _h[1]->fill(recoilW(p, PID::D0BAR));
      	}
      	else if(p.pid()==PID::B0) {
      	  _nB[1]->fill();
      	  if (isSemileptonicDecay(p, {PID::DMINUS,PID::POSITRON,PID::NU_E}))  _h[2]->fill(recoilW(p, PID::DMINUS));
      	  if (isSemileptonicDecay(p, {PID::DMINUS,PID::ANTIMUON,PID::NU_MU})) _h[3]->fill(recoilW(p, PID::DMINUS));
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double tau[2] = {1638e-15,1519e-15};
      double hbar = 6.582119569e-10;
      for(unsigned int ix=0;ix<4;++ix) {
      	double fact = ix<2 ? hbar/tau[0] : hbar/tau[1];
      	CounterPtr c = ix<2 ? _nB[0] : _nB[1];
      	scale(_h[ix],fact/ *c);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    CounterPtr _nB[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2169621);

}
