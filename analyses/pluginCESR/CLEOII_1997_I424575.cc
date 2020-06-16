// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Sigma*_c spectra
  class CLEOII_1997_I424575 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1997_I424575);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_cTheta,2,1,1);
      book(_h_x     ,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.518));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==4114 or Cuts::abspid==4224)) {
	// momentum fraction
	double xp = p.momentum().p3().mod()/Pmax;
        _h_x->fill(xp);
	if(p.children().size()!=2) continue;
	// decay angle
	Particle pi;
	int sign = p.pid()/p.abspid();
	if(p.abspid()==4224) {
	  if(p.children()[0].pid() == sign*4122 &&
	     p.children()[1].pid() == sign*211 ) {
	    pi=p.children()[1];
	  }
	  else if(p.children()[1].pid() == sign*4122 &&
		  p.children()[0].pid() == sign*211 ) {
	    pi=p.children()[0];
	  }
	  else
	    continue;
	}
	else {
	  if(p.children()[0].pid() == sign*4122 &&
	     p.children()[1].pid() == -sign*211 ) {
	    pi=p.children()[1];
	  }
	  else if(p.children()[1].pid() == sign*4122 &&
		  p.children()[0].pid() == -sign*211 ) {
	    pi=p.children()[0];
	  }
	  else
	    continue;
	}
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axis = boost.transform(pi.momentum()).p3().unit();
	double cosL  = axis.dot(p.momentum().p3().unit());
	_h_cTheta->fill(cosL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x     ,1.,false);
      normalize(_h_cTheta,1.,false);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_h_cTheta;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1997_I424575);

}
