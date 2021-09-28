// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {



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
    if (!ws.size())  return std::valarray<double>{1.0};
    size_t N = _weightIndices.size();
    if (N == ws.size())  return ws;
    std::valarray<double> new_ws(N);
    for (size_t i = 0; i < N; ++i) {
      new_ws[i] = ws[_weightIndices[i]];
    }
    return new_ws;
  }


}
