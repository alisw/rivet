// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectrum for Omega_c
  class BABAR_2007_I746745 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2007_I746745);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(),"UFS");
      book(_h_p,1,1,1);
      book(_b  ,2,1,1);
      book(_r  ,3,1,1);
      book(_ups,"/TMP/ups");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idOmega = 4332;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      bool ups = !ufs.particles(Cuts::pid==300553).empty();
      if(ups) _ups->fill();
      for (const Particle& p : ufs.particles(Cuts::abspid==idOmega)) {
	_h_p->fill(p.momentum().p3().mod());
	if(p.children().size()==2) {
	  int sign = p.pid()/p.abspid();
	  if((p.children()[0].pid()==sign*3334 &&
	      p.children()[1].pid()==sign*211) ||
	     (p.children()[1].pid()==sign*3334 &&
	      p.children()[0].pid()==sign*211) ) {
	    if(ups) _b->fill(0.5);
	    else    _r->fill(0.5);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_p);
      if(_ups->effNumEntries()!=0) {
	scale(_b,0.5/ *_ups);
      }
      scale(_r,crossSection()/sumOfWeights()/femtobarn);
    }
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_p,_b,_r;
    CounterPtr _ups;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2007_I746745);

}
