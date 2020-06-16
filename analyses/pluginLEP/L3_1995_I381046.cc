// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B* production
  class L3_1995_I381046 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_1995_I381046);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // // Initialise and register projections
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_bStar, "/TMP/cbStar ");
      book(_c_B    , "/TMP/cB     ");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==513 or Cuts::abspid==523 or
									       Cuts::abspid==511 or Cuts::abspid==521)) {
	// count number of Bs not from mixing or B*
	if(p.abspid()==511 || p.abspid()==521) {
	  if(p.parents()[0].abspid()==p.abspid()) continue;
	  if(p.parents()[0].abspid()==513 || p.parents()[0].abspid()==523) continue;
	  _c_B->fill(); 
	}
	// B*
	else {
	  _c_bStar->fill();
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // no of B*/B+B*
      Scatter2DPtr h1;
      book(h1,1,1,1);
      Counter ctemp = *_c_bStar+*_c_B;
      double val = _c_bStar->val()/ctemp.val();
      double err = val*sqrt(sqr(_c_bStar->err()/_c_bStar->val())+sqr(ctemp.err()/ctemp.val()));
      h1->addPoint(91.2,val,make_pair(0.5,0.5),make_pair(err,err) );
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_bStar,_c_B;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_1995_I381046);


}
