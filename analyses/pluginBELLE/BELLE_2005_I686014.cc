// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief charm hadron production
  class BELLE_2005_I686014 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2005_I686014);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histos
      if(fuzzyEquals(sqrtS(),10.52,1e-3))
	_mode=1;
      else if(fuzzyEquals(sqrtS(),10.58,1e-3))
	_mode=2;
      else
	MSG_ERROR("Beam energy not supported!");
      for(unsigned int ix=0;ix<7;++ix) {
	if(_mode==1)
	  book(_r[ix],2,1,ix+1);
	else
	  book(_r[ix],1,1,ix+1);
	book(_h[ix],2+_mode,1,ix+1);
      }
      book(_c,"TMP/wgt");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // unstable particles
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      if(_mode==2 && ufs.particles(Cuts::pid==300553).size()!=1)
	vetoEvent;
      _c->fill();
      for(const Particle & p : ufs.particles(Cuts::abspid==411 ||
					     Cuts::abspid==421 ||
					     Cuts::abspid==431 ||
					     Cuts::abspid==413 ||
					     Cuts::abspid==423 ||
					     Cuts::abspid==4122 )) {
	double pmax = sqrt(0.25*sqr(sqrtS())-sqr(p.mass()));
	double xp = p.momentum().p3().mod()/pmax;
	if(p.abspid()==421) {
	  _r[0]->fill(0.5);
	  _h[0]->fill(xp);
	}
	else if(p.abspid()==421) {
	  _r[0]->fill(0.5);
	  _h[0]->fill(xp);
	}
	else if(p.abspid()==411) {
	  _r[1]->fill(0.5);
	  _h[1]->fill(xp);
	}
	else if(p.abspid()==431) {
	  _r[2]->fill(0.5);
	  _h[2]->fill(xp);
	}
	else if(p.abspid()==4122) {
	  _r[3]->fill(0.5);
	  _h[3]->fill(xp);
	}
	else if(p.abspid()==413) {
	  _r[4]->fill(0.5);
	  _h[4]->fill(xp);
	  _r[5]->fill(0.5);
	  _h[5]->fill(xp);
	}
	else if(p.abspid()==423) {
	  _r[6]->fill(0.5);
	  _h[6]->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_mode==1) {
	for(unsigned int ix=0;ix<7;++ix) {
	  scale(_r[ix],crossSection()/picobarn/sumOfWeights());
	  scale(_h[ix],crossSection()/nanobarn/sumOfWeights());
	}
      }
      else {
	for(unsigned int ix=0;ix<7;++ix) {
	  scale(_r[ix], 0.5/ *_c);
	  scale(_h[ix],crossSection()/nanobarn/sumOfWeights());
	}
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h[7];
    Histo1DPtr _r[7];
    CounterPtr _c;
    unsigned int _mode=0;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2005_I686014);

}
