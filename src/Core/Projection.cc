// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Cmp.hh"

namespace Rivet {


  Projection::Projection()
    : _name("BaseProjection"), _isValid(true)
  {
    addPdgIdPair(PID::ANY, PID::ANY);
  }


  Projection:: ~Projection() = default;


  Projection& Projection::operator = (const Projection&) { return *this; }


  bool Projection::before(const Projection& p) const {
    // actual ordering is irrelevant, return true if projections not equal
    const std::type_info& thisid = typeid(*this);
    const std::type_info& otherid = typeid(p);
    if (thisid == otherid) {
      const CmpState cmpst = compare(p);
      const bool cmp = (cmpst != CmpState::EQ);
      MSG_TRACE("Comparing projections of same RTTI type: " << this << " < " << &p << " = " << cmp);
      return cmp;
    } else {
      const bool cmp = thisid.before(otherid);
      MSG_TRACE("Ordering projections of different RTTI type: " << this << " < " << &p << " = " << cmp);
      return cmp;
    }
  }


  const set<PdgIdPair> Projection::beamPairs() const {
    set<PdgIdPair> ret = _beamPairs;
    set<ConstProjectionPtr> projs = getProjections();
    for (set<ConstProjectionPtr>::const_iterator ip = projs.begin(); ip != projs.end(); ++ip) {
      ConstProjectionPtr p = *ip;
      getLog() << Log::TRACE << "Proj addr = " << p << '\n';
      if (p) ret = intersection(ret, p->beamPairs());
    }
    return ret;
  }


  Cmp<Projection> Projection::mkNamedPCmp(const Projection& otherparent, const string& pname) const {
    return pcmp(*this, otherparent, pname);
  }

  Cmp<Projection> Projection::mkPCmp(const Projection& otherparent, const string& pname) const {
    return pcmp(*this, otherparent, pname);
  }


}
