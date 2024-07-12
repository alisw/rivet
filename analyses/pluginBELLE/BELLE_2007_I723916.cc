// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c(2880) -> Sigma_c pi
  class BELLE_2007_I723916 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I723916);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      // set the PDG code
      _pid = getOption<double>("PID", 4126);
      // projections
      declare(UnstableParticles(Cuts::abspid==_pid), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // Get beams and average beam momentum
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	// xp cut
	double xp = p.momentum().p3().mod()/sqrt(0.25*sqrtS()-p.mass2());
	if (xp<.7) continue;
	// find the decay products
	int sign = p.pid()/p.abspid();
	// first decay
	if(p.children().size()!=2) continue;
	Particle Sigma,pi1;
	if((p.children()[0].abspid()== 4222*sign &&
	    p.children()[1].abspid()==-211 *sign) ||
	   (p.children()[0].abspid()== 4112*sign &&
	    p.children()[1].abspid()== 211 *sign) ) {
	  Sigma = p.children()[0];
	  pi1   = p.children()[1];
	}
	else if((p.children()[1].abspid()== 4222*sign &&
		 p.children()[0].abspid()==-211 *sign) ||
		(p.children()[1].abspid()== 4112*sign &&
		 p.children()[0].abspid()== 211 *sign) ) {
	  Sigma = p.children()[1];
	  pi1   = p.children()[0];
	}
	else
	  continue;
	Vector3 axis1 = p.momentum().p3().unit();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axis2 = boost.transform(pi1.momentum()).p3().unit();
	_h[0]->fill(axis1.dot(axis2));
	Vector3 axis3 = axis.cross(axis1).unit();
	Vector3 axis4 = axis1.cross(axis2).unit();
	double phi = acos(axis3.dot(axis4));
	_h[1]->fill(phi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I723916);

}
