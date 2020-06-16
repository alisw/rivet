// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /*
  double Event::weight() const {
    // Get the weight index, which defaults to 0, i.e. nominal
    // NB. This should normally only perform the slow env var lookup once per run
    // NB. An explicitly set index of -1 is used to ignore all weights, for debugging
    static int WEIGHT_INDEX = -2;
    if (WEIGHT_INDEX < -1) {
      WEIGHT_INDEX = getEnvParam<size_t>("RIVET_WEIGHT_INDEX", 0);
      Log::getLog("Core.Weight") << Log::TRACE << "Got weight index from env/default = "<< WEIGHT_INDEX << endl;
    }
    // If RIVET_WEIGHT_INDEX=-1, or there are no event weights, return 1
    if (WEIGHT_INDEX == -1 || genEvent()->weights().empty()) return 1.0;
    // Otherwise return the appropriate weight index
    return _genevent.weights()[WEIGHT_INDEX];
  }
  */


  ParticlePair Event::beams() const { return Rivet::beams(*this); }


  double Event::sqrtS() const { return Rivet::sqrtS(beams()); }


  double Event::asqrtS() const { return Rivet::asqrtS(beams()); }


  void Event::_init(const GenEvent& ge) {
    // Use Rivet's preferred units if possible
    #ifdef HEPMC_HAS_UNITS
    _genevent.use_units(HepMC::Units::GEV, HepMC::Units::MM);
    #endif
  }


  void Event::_strip(GenEvent & ge) {
    HepMCUtils::strip(ge);
  }


  const Particles& Event::allParticles() const {
    if (_particles.empty()) { //< assume that empty means no attempt yet made
      for (ConstGenParticlePtr gp : HepMCUtils::particles(genEvent())) {
        _particles += Particle(gp);
      }
    }
    return _particles;
  }


  std::valarray<double> Event::weights() const {
    const std::valarray<double> ws = HepMCUtils::weights(_genevent);
    return ws.size() > 0 ? ws : std::valarray<double>{1.0};
  }


}
