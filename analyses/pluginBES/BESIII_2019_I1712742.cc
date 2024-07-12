// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s -> eta, eta' semi-leptonic
  class BESIII_2019_I1712742 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1712742);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],ix+1,1,iy+1);
      // counter for normalization
      book(_nDs, "TMP/nDs");
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
      // Loop over Ds mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==431)) {
	_nDs->fill();
        if (isSemileptonicDecay(p, {PID::ETA, PID::POSITRON, PID::NU_E})) {
          _h[0][0]->fill(q2(p, PID::ETA));
          _h[0][1]->fill(q2(p, PID::ETA));
        }
	else if(isSemileptonicDecay(p, {PID::ETAPRIME, PID::POSITRON, PID::NU_E})) {
          _h[1][0] ->fill(q2(p, PID::ETAPRIME));
          _h[1][1] ->fill(q2(p, PID::ETAPRIME));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // D_s lifetime pdg2018 (ns)
      double tau = 504e-6;
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  scale(_h[ix][iy],1./tau/ *_nDs);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nDs;
    Histo1DPtr _h[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1712742);

}
