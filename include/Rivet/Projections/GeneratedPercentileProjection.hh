// -*- C++ -*-
#ifndef RIVET_GENERATEDPERCENTILEPROJECTION_HH
#define RIVET_GENERATEDPERCENTILEPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {

class GeneratedPercentileProjection: public SingleValueProjection {
public:
  
  GeneratedPercentileProjection() {
    setName("GeneratedPercentileProjection");
    declare(HepMCHeavyIon(), "HepMC");
  }

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(GeneratedPercentileProjection);

protected:

  void project(const Event& e) {
    clear();
    set(apply<HepMCHeavyIon>(e, "HepMC").centrality());
   }

  CmpState compare(const Projection& p) const {
    return CmpState::EQ;
  }
  
};

}

#endif
