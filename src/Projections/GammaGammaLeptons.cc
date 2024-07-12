// -*- C++ -*-
#include "Rivet/Projections/GammaGammaLeptons.hh"

namespace Rivet {


  CmpState GammaGammaLeptons::compare(const Projection& p) const {
    const GammaGammaLeptons& other = pcast<GammaGammaLeptons>(p);
    return mkNamedPCmp(other, "Beam") || mkNamedPCmp(other, "LFS") ||
      mkNamedPCmp(other, "IFS") || cmp(_sort, other._sort);
  }


  void GammaGammaLeptons::project(const Event& e) {
    // Find incoming lepton beams
    _incoming = applyProjection<Beam>(e, "Beam").beams();
    // need two leptonic beams
    if(! PID::isLepton(_incoming. first.pid()) ||
       ! PID::isLepton(_incoming.second.pid()) ) {
      fail();
      return;
    }

    // If no graph-connected scattered lepton, use the hardest
    // (preferably same-flavour) prompt FS lepton in the event.
    const FinalState & fs = applyProjection<FinalState>(e, "LFS");
    Particles fsleptons;
    if ( _sort == ET )
       fsleptons = fs.particles(isLepton, cmpMomByEt);
    else if ( _sort == ETA && _incoming.first.momentum().pz() >= 0.0 )
      fsleptons = fs.particles(isLepton, cmpMomByDescEta);
    else if ( _sort == ETA && _incoming.first.momentum().pz() < 0.0 )
      fsleptons = fs.particles(isLepton, cmpMomByEta);
    else
      fsleptons = fs.particles(isLepton, cmpMomByE);

    // loop over the two beam particles
    for(unsigned int ix=0;ix<2;++ix) {
      Particle inc = ix == 0 ? _incoming.first : _incoming.second;
      // resort if required
      if(ix==1 && _sort ==ETA ) {
	if ( _sort == ETA && _incoming.second.momentum().pz() >= 0.0 )
	  sort(fsleptons.begin(),fsleptons.end(), cmpMomByDescEta);
	else if ( _sort == ETA && _incoming.second.momentum().pz() < 0.0 )
	  sort(fsleptons.begin(),fsleptons.end(), cmpMomByEta);
      }
      Particles sfleptons =
	filter_select(fsleptons, Cuts::pid == inc.pid());
      if ( sfleptons.empty() ) sfleptons = fsleptons;

      if ( _isolDR > 0.0 ) {
	const Particles & other =
	  applyProjection<FinalState>(e, "IFS").particles();
	while (!sfleptons.empty()) {
	  bool skip = false;
	  Particle testlepton = sfleptons.front();
	  for ( auto p: other ) {
	    if ( skip ) break;
	    if ( deltaR(p, testlepton) < _isolDR ) skip = true;
	    for ( auto c : testlepton.constituents() ) {
	      if ( c.genParticle() == p.genParticle() ) {
		skip = false;
		break;
	      }
	    }
	  }
	  if ( !skip ) break;
	  sfleptons.erase(sfleptons.begin());
	}
      }
      if ( !sfleptons.empty() ) {
	if(ix==0) {
	  _outgoing.first = sfleptons.front();
	}
	else {
	  _outgoing.second = sfleptons.front();
	}
      }
      else {
	fail();
      }
      
    }
  }


}
