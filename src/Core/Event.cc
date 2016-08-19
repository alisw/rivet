// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  ParticlePair Event::beams() const { return Rivet::beams(*this); }

  // PdgIdPair Event::beamIds() const { return pids(beams()); }

  double Event::sqrtS() const { return Rivet::sqrtS(beams()); }

  double Event::asqrtS() const { return Rivet::asqrtS(beams()); }

  // Vector3 Event::beamCMSBoost() const { return Rivet::beamCMSBoost(*this); }

  // LorentzTransform Event::beamCMSTransform() const { return Rivet::beamCMSTransform(*this); }



  void Event::_init(const GenEvent& ge) {
    // Use Rivet's preferred units if possible
    #ifdef HEPMC_HAS_UNITS
    _genevent.use_units(HepMC::Units::GEV, HepMC::Units::MM);
    #endif

    // Use the conventional alignment
    // _geNormAlignment();

    /// @todo Filter the event to remove generator-specific particles: optional
    /// behaviour? Maybe disableable in an inconvenient way, e.g. with an env
    /// var, to communicate the appropriate distaste for this sort of truth
    /// analysis ;-)

    // Debug printout to check that copying/mangling has worked
    /// @todo Enable this when HepMC has been fixed to allow printing to a stream like the Rivet logger.
    //_genevent.print();
  }


  // namespace { // unnamed namespace for hiding
  //
  //   void _geRot180x(GenEvent& ge) {
  //     /// @todo Use nicer iterators over HepMC particles
  //     for (HepMC::GenEvent::particle_iterator ip = ge.particles_begin(); ip != ge.particles_end(); ++ip) {
  //       const HepMC::FourVector& mom = (*ip)->momentum();
  //       (*ip)->set_momentum(HepMC::FourVector(mom.px(), -mom.py(), -mom.pz(), mom.e()));
  //     }
  //     /// @todo Use nicer iterators over HepMC vertices
  //     for (HepMC::GenEvent::vertex_iterator iv = ge.vertices_begin(); iv != ge.vertices_end(); ++iv) {
  //       const HepMC::FourVector& pos = (*iv)->position();
  //       (*iv)->set_position(HepMC::FourVector(pos.x(), -pos.y(), -pos.z(), pos.t()));
  //     }
  //   }
  //
  // }


  // void Event::_geNormAlignment() {
  //   if (!_genevent.valid_beam_particles()) return;
  //   typedef pair<HepMC::GenParticle*, HepMC::GenParticle*> GPPair;
  //   GPPair bps = _genevent.beam_particles();
  //
  //   // Rotate e+- p and ppbar to put p along +z
  //   /// @todo Is there an e+ e- convention for longitudinal boosting, e.g. at B-factories? Different from LEP?
  //   // if (compatible(beamids, make_pdgid_pair(ELECTRON, PROTON)) ||
  //   //     compatible(beamids, make_pdgid_pair(POSITRON, PROTON)) ||
  //   //     compatible(beamids, make_pdgid_pair(ANTIPROTON, PROTON)) ) {
  //   //   Log::getLog("Rivet.Event") << Log::TRACE << "May need to rotate event..." << endl;
  //   bool rot = false;
  //   const HepMC::GenParticle* plusgp = 0;
  //   if (bps.first->pdg_id() != PID::PROTON || bps.second->pdg_id() != PID::PROTON) {
  //     if (bps.first->pdg_id() == PID::PROTON) {
  //       plusgp = bps.first;
  //     } else if (bps.second->pdg_id() == PID::PROTON) {
  //       plusgp = bps.second;
  //     }
  //     if (plusgp && plusgp->momentum().pz() < 0) {
  //       rot = true;
  //     }
  //   }
  //
  //   // Do the rotation
  //   if (rot) {
  //     if (Log::getLog("Rivet.Event").isActive(Log::TRACE)) {
  //       Log::getLog("Rivet.Event") << Log::TRACE << "Rotating event\n"
  //                                  << "Before rotation: "
  //                                  << bps.first->pdg_id() << "@pz=" << bps.first->momentum().pz()/GeV << ", "
  //                                  << bps.second->pdg_id() << "@pz=" << bps.second->momentum().pz()/GeV << endl;
  //     }
  //     _geRot180x(_genevent);
  //   }
  // }


  const Particles& Event::allParticles() const {
    if (_particles.empty()) { //< assume that empty means no attempt yet made
      for (const GenParticle* gp : particles(genEvent())) {
        _particles += Particle(gp);
      }
    }
    return _particles;
  }


  double Event::weight() const {
    return (!_genevent.weights().empty()) ? _genevent.weights()[0] : 1.0;
  }


}
