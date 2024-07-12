// -*- C++ -*-
#include "Rivet/Projections/DISDiffHadron.hh"

namespace Rivet {


  CmpState DISDiffHadron::compare(const Projection& p) const {
    return mkNamedPCmp(p, "Beam") || mkNamedPCmp(p, "FS");
  }


  void DISDiffHadron::project(const Event& e) {

      // Find incoming lepton beam
      const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();
      bool firstIsHadron = PID::isHadron(inc.first.pid());
      bool secondIsHadron = PID::isHadron(inc.second.pid());
      if (firstIsHadron && !secondIsHadron) {
        _incoming = inc.first;
      } else if (!firstIsHadron && secondIsHadron) {
        _incoming = inc.second;
      } else {
        fail();
        return;
      }

      const FinalState & fs = applyProjection<FinalState>(e, "FS");
      Particles fshadrons;
      if ( _incoming.momentum().pz() >= 0.0 )
        fshadrons = fs.particles(isHadron, cmpMomByDescEta);
      else
        fshadrons = fs.particles(isHadron, cmpMomByEta);

      Particles sfhadrons = filter_select(fshadrons,
                                          Cuts::pid == _incoming.pid());
      MSG_DEBUG("SF hadrons = " << sfhadrons.size() <<
                ", all hadrons = " << fshadrons.size());
      if (!sfhadrons.empty()) {
        _outgoing = sfhadrons.front();
      } else if (!fshadrons.empty()) {
        _outgoing = fshadrons.front();
      } else {
        fail();
      }

    }


}
