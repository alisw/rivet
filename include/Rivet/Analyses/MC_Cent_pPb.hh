// -*- C++ -*-
#ifndef RIVET_MC_Cent_pPb_HH
#define RIVET_MC_Cent_pPb_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/TriggerProjection.hh"

namespace Rivet {

/// Example of a centrality observable projection for pPb that uses
/// summed Et in the Pb direction.
class MC_SumETFwdPbCentrality: public SingleValueProjection {

public:

  /// Constructor.
  MC_SumETFwdPbCentrality() {
    declare(FinalState(Cuts::eta < -3.2 && Cuts::eta > -4.9 && Cuts::pT > 0.1*GeV),
	    "FSSumETFwdCentrality");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(MC_SumETFwdPbCentrality);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    clear();
    const FinalState & fsfwd =
      apply<FinalState>(e, "FSSumETFwdCentrality");
    double estimate = 0.0;
    for ( const Particle & p : fsfwd.particles() ) {
      estimate += p.Et();
    }
    set(estimate);
  }
  
  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "FSSumETFwdCentrality");
  }

};
    
/// Example of a trigger projection for minimum bias pPb requiring at
/// least one charged particle in both forward and backward direction.
class MC_pPbMinBiasTrigger: public TriggerProjection {

public:

  /// Constructor.
  MC_pPbMinBiasTrigger() {
    declare(FinalState(Cuts::eta < -3.2 && Cuts::eta > -4.9 && Cuts::pT > 0.1*GeV),
	    "FSSumETFwdCentrality");
      declare(ChargedFinalState(Cuts::eta > 2.09 && Cuts::eta < 3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBB");
      declare(ChargedFinalState(Cuts::eta < -2.09 && Cuts::eta > -3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBF");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(MC_pPbMinBiasTrigger);

protected:

  /// Perform the projection on the Event
  void project(const Event& event) {
    pass();
      if ( applyProjection<FinalState>(event,"MBF").particles().empty() ||
	   applyProjection<FinalState>(event,"MBB").particles().empty() )
        fail();
  }
  
  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "MBF") || mkNamedPCmp(p, "MBB");
  }

};
    

}

#endif
