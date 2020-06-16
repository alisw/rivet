// -*- C++ -*-
#ifndef RIVET_ATLAS_COMMON_HH
#define RIVET_ATLAS_COMMON_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/TriggerProjection.hh"

namespace Rivet {
/// Common projections for ATLAS trigger conditions and centrality.

  namespace ATLAS {


    /// Centrality projection for pPb collisions (one sided)
class SumET_PB_Centrality: public SingleValueProjection {

public:

  /// Constructor.
    SumET_PB_Centrality() {
    declare(FinalState(Cuts::eta < -3.2 && Cuts::eta > -4.9 && Cuts::pT > 0.1*GeV),
	    "SumET_PB_Centrality");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SumET_PB_Centrality);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    clear();
    const FinalState & fsfwd =
      apply<FinalState>(e, "SumET_PB_Centrality");
    double estimate = 0.0;
    for ( const Particle & p : fsfwd.particles() ) {
      estimate += p.Et();
    }
    set(estimate);
  }

  /// Compare projections
  CmpState compare(const Projection& p) const {
    return mkNamedPCmp(p, "SumET_PB_Centrality");
  }

};


/// Centrality projection for PbPb collisions (two sided)
class SumET_PBPB_Centrality: public SingleValueProjection {

public:

  /// Constructor.
    SumET_PBPB_Centrality() {
    declare(FinalState(Cuts::abseta > 3.2 && Cuts::abseta < 4.9 && Cuts::pT > 0.1*GeV),
	    "SumET_PBPB_Centrality");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SumET_PBPB_Centrality);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    clear();
    const FinalState & fsfwd =
      apply<FinalState>(e, "SumET_PBPB_Centrality");
    double estimate = 0.0;
    for ( const Particle & p : fsfwd.particles() ) {
      estimate += p.Et();
    }
    set(estimate);
  }

  /// Compare projections
  CmpState compare(const Projection& p) const {
    return mkNamedPCmp(p, "SumET_PBPB_Centrality");
  }

};

/// ATLAS min bias trigger conditions.
class MinBiasTrigger: public TriggerProjection {

public:

  /// Constructor.
  MinBiasTrigger() {
      declare(ChargedFinalState(Cuts::eta > 2.09 && Cuts::eta < 3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBB");
      declare(ChargedFinalState(Cuts::eta < -2.09 && Cuts::eta > -3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBF");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(MinBiasTrigger);

protected:

  /// Perform the projection on the Event
  void project(const Event& event) {
    pass();
      if ( applyProjection<FinalState>(event,"MBF").particles().empty() ||
	   applyProjection<FinalState>(event,"MBB").particles().empty() )
        fail();
  }

  /// Compare projections
  CmpState compare(const Projection& p) const {
    return mkNamedPCmp(p, "MBF") || mkNamedPCmp(p, "MBB");
  }

};


}
}

#endif
