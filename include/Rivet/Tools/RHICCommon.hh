// -*- C++ -*-
#ifndef RIVET_ATLAS_COMMON_HH
#define RIVET_ATLAS_COMMON_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/TriggerProjection.hh"

namespace Rivet {
/// Common projections for RHIC experiments' trigger conditions and centrality.

class STAR_BES_Centrality : public SingleValueProjection {
public:
  STAR_BES_Centrality() {
    declare(ChargedFinalState(Cuts::abseta < 0.5 &&
      Cuts::absrap < 0.1 && Cuts::pT > 0.2 * GeV),
      "STAR_BES_Centrality");
  }
  
  // Destructor
  virtual ~STAR_BES_Centrality() {}
  
  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(STAR_BES_Centrality);

protected:
  void project(const Event& e) {
    clear();
    double estimate = 
      apply<FinalState>(e, "STAR_BES_Centrality").particles().size();
		set(estimate);
  }

  /// Compare projections
  virtual CmpState compare(const Projection& p) const {
    return mkNamedPCmp(p, "STAR_BES_Centrality");
  }
};


/// @brief BRAHMS Centrality projection.
class BRAHMSCentrality : public SingleValueProjection {
public:
  // Constructor
  BRAHMSCentrality() : SingleValueProjection() {
    // Using here the BRAHMS reaction centrality from eg. 1602.01183, which
    // might not be correct.
    declare(ChargedFinalState(Cuts::pT > 0.1*GeV && Cuts::abseta < 2.2),
      "ChargedFinalState");
  }
  // Destructor
  virtual ~BRAHMSCentrality() {}

  // Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(BRAHMSCentrality);

protected:
  // Do the projection. Count the number of charged particles in
  // the specified range.
  virtual void project(const Event& e) {
    clear();
    set(apply<ChargedFinalState>
      (e, "ChargedFinalState").particles().size());
  }

  // Compare to another projection.
  virtual CmpState compare(const Projection& p) const {
    return mkNamedPCmp(p, "BRAHMSCentrality");
  }

};
}
#endif
