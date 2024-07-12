// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B_s0 -> D_s* semileptonic decay
  class LHCB_2020_I1787090 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2020_I1787090);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_w , 1, 1, 1);
    }

    // Calculate the Q2 using mother and daugher meson
    double w(const Particle& B, int mesonID) {
      Particle D = filter_select(B.children(), Cuts::pid==mesonID)[0];
      FourMomentum q = B.mom() -D.mom() ;
      double q2 = q*q;
      return 0.5*(sqr(B.mass())+sqr(D.mass())-q2)/B.mass()/D.mass();
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
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==531)) {
	if(p.pid()<0) {
	  if(isSemileptonicDecay(p, {433, PID::MUON, PID::NU_MUBAR})) {
	    _h_w->fill(w(p, 433));
	  }
	}
	else {
	  if(isSemileptonicDecay(p, {-433, PID::ANTIMUON, PID::NU_MU})) {
	    _h_w->fill(w(p, -433));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize to unity
      normalize(_h_w);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_w;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2020_I1787090);

}
