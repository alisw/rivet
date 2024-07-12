// -*- C++ -*-
#include "Rivet/Projections/GammaGammaFinalState.hh"

namespace Rivet {


  void GammaGammaFinalState::project(const Event& e) {
    const GammaGammaKinematics& ggkin = apply<GammaGammaKinematics>(e, "Kinematics");
    if ( ggkin.failed() ) {
      fail();
      return;
    }

    const GammaGammaLeptons& gglep = ggkin.apply<GammaGammaLeptons>(e, "Lepton");
    if ( ggkin.failed() ) {
      fail();
      return;
    }

    const FinalState& fs = apply<FinalState>(e, "FS");

    // Fill the particle list with all particles _other_ than the GammaGamma scattered
    // lepton, with momenta boosted into the appropriate frame.
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size()-1);
    ConstGenParticlePtr lep1 = gglep.out().first .genParticle();
    ConstGenParticlePtr lep2 = gglep.out().second.genParticle();
    // Ensure that we skip the GammaGamma leptons
    for (const Particle& p : fs.particles()) { 
      if (p.genParticle() != lep1 &&
	  p.genParticle() != lep2)  _theParticles.push_back(p);
    }
  }


}
