// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// @brief D*+/- polarization at 29 GeV
  class HRS_1987_I250823 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1987_I250823);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(Sphericity(FinalState()), "Sphericity");
      declare(UnstableParticles(), "UFS" );
      for(unsigned int i=0;i<9;++i) {
       	unsigned int ix(0),iy(1);
       	if(i<3) {
       	  ix=1;
       	  iy=i+1;
       	}
       	else if(i==3) {
       	  ix=2;
       	}
       	else if(i<6) {
       	  ix=3;
       	  iy=i-3;
       	}
       	else {
       	  ix=i-2;
       	}
	book(_p_rho00[i],ix,iy,1);
	book(_p_rho11[i],ix,iy,2);
	book(_p_rho10[i],ix,iy,3);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      
      // sphericity, to define an axis
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      // loop over the particles
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==413)) {
	int sign = p.pid()/p.abspid();
	Particle decay;
	double xE = p.momentum().E()/meanBeamMom;
	if(p.children()[0].pid()==sign*421 && 
	   p.children()[1].pid()==sign*211) {
	  decay = p.children()[1];
	}
	else if(p.children()[1].pid()==sign*421 && 
		p.children()[0].pid()==sign*211) {
	  decay = p.children()[0];
	}
	else
	  continue;
	// axis and ctheta
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 e1z = p.p3().unit();	
	FourMomentum pp = boost.transform(decay.momentum());
	Vector3 axis1 = boost.transform(decay.momentum()).p3().unit();
	double ctheta = e1z.dot(axis1);
	// y and z axis
	Vector3 e1y = e1z.cross(axis).unit();
	Vector3 e1x = e1y.cross(e1z).unit();
	double phi = atan2(e1y.dot(axis1),e1x.dot(axis1));
	double w1 = 0.5*(5.*sqr(ctheta)-1.);
	double w2 = -1.25*(1.-sqr(ctheta))*cos(2.*phi);
	double w3 = -1.25*sqrt(2)*ctheta*sqrt(1.-sqr(ctheta))*cos(phi);
	// fill the hists by x_E
	for(unsigned int ix=0;ix<6;++ix) {
	  _p_rho00[ix]->fill(xE,w1);
	  _p_rho11[ix]->fill(xE,w2);
	  _p_rho10[ix]->fill(xE,w3);
	}
	_p_rho00[6]->fill(1.,w1);
	_p_rho11[6]->fill(1.,w2);
	_p_rho10[6]->fill(1.,w3);
	// using jet axis
	double pT = sqrt(sqr(sphericity.sphericityMajorAxis().dot(p.momentum().p3()))+
			 sqr(sphericity.sphericityMinorAxis().dot(p.momentum().p3())));
	Vector3 axis2;
	if(p.momentum().p3().dot(sphericity.sphericityAxis())>=0.) {
	  axis2 = sphericity.sphericityAxis();
	}
	else {
	  axis2 =-sphericity.sphericityAxis();
	}
	Vector3 e2y = e1z.cross(axis2).unit();
	Vector3 e2x = e2y.cross(e1z).unit();
	// alpha and beta
	phi = atan2(e2y.dot(axis1),e2x.dot(axis1));
	w1 = 0.5*(5.*sqr(ctheta)-1.);
	w2 = -1.25*(1.-sqr(ctheta))*cos(2.*phi);
	w3 = -1.25*sqrt(2)*ctheta*sqrt(1.-sqr(ctheta))*cos(phi);
	if(xE<0.4) continue;
	_p_rho00[7]->fill(1.,w1);
	_p_rho11[7]->fill(1.,w2);
	_p_rho10[7]->fill(1.,w3);
	if(pT<0.75) {
	  _p_rho00[7]->fill(2.,w1);
	  _p_rho11[7]->fill(2.,w2);
	  _p_rho10[7]->fill(2.,w3);
	}
	else {
	  _p_rho00[7]->fill(3.,w1);
	  _p_rho11[7]->fill(3.,w2);
	  _p_rho10[7]->fill(3.,w3);
	}
	if(sphericity.sphericity()>=0.1) {
	  _p_rho00[8]->fill(1.,w1);
	  _p_rho11[8]->fill(1.,w2);
	  _p_rho10[8]->fill(1.,w3);
	}
	else {
	  _p_rho00[8]->fill(2.,w1);
	  _p_rho11[8]->fill(2.,w2);
	  _p_rho10[8]->fill(2.,w3);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    //@}


    /// @name Histograms
    //@{
    Profile1DPtr _p_rho00[9],_p_rho11[9],_p_rho10[9];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(HRS_1987_I250823);


}
