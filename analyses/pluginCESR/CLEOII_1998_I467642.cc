// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Upsilon(2S) -> Upsilon(1S) pi+pi-
  class CLEOII_1998_I467642 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1998_I467642);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      book(_h_charged_ex   ,1,1,1);
      book(_h_charged_inc  ,2,1,1);
      book(_h_neutral_ex   ,3,1,1);
      book( _h_lept_ctheta ,4,1,1);
      book(_h_lept_phi     ,4,1,2);
      book(_h_piStar_ctheta,5,1,1);
      book(_h_piStar_phi   ,5,1,2);
      book(_h_pi_ctheta    ,6,1,1);
      book(_h_pi_phi       ,6,1,2);

    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles& pi0, Particles & onium) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
       	else if (id == PID::PI0) {
       	  pi0.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)%1000==443 || abs(id)%1000==553) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,onium);
	}
	else
	  ++nstable;
      }
    }

    void findLeptons(const Particle & mother,
		     unsigned int & nstable,
		     Particles& lp, Particles& lm) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == 11 || id == 13 ) {
	  lm.push_back(p);
	  ++nstable;
	}
       	else if (id == -11 || id==-13) {
       	  lp.push_back(p);
       	  ++nstable;
       	}
	else if ( !p.children().empty() ) {
	  findLeptons(p,nstable,lp,lm);
	}
	else
	  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==100553)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, onium;
	findDecayProducts(ups,nstable,pip,pim,pi0,onium);
	// check for onium
	if(onium.size() !=1 || onium[0].pid()!=553 || nstable !=3) continue;
	//pi+pi-
	if( pip.size()==1 && pim.size() ==1) {
	  FourMomentum q = pip[0].momentum()+pim[0].momentum();
	  _h_charged_ex->fill(q.mass()/GeV);
	  _h_charged_inc->fill(q.mass()/GeV);
	  // leptons from Upsilon(1S) decay
	  nstable = 0;
	  Particles lp, lm;
	  findLeptons(onium[0],nstable,lp,lm);
	  if(nstable==2&&lp.size()==1&&lm.size()==1) {
	    _h_lept_ctheta->fill(cos(lp[0].momentum().polarAngle()));
	    _h_lept_phi   ->fill(lp[0].momentum().azimuthalAngle());
	  }
	  // pions in rest frame
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(q.betaVec());
	  FourMomentum ppi = boost.transform(pip[0].momentum());
	  Vector3 axis1 = q.p3().unit();
	  Vector3 axis2 = Vector3(0.,0.,1.).cross(axis1).unit();
	  Vector3 axis3 = axis1.cross(axis2).unit();
	  _h_piStar_ctheta->fill(axis1.dot(ppi.p3().unit()));
	  double phi = atan2(axis3.dot(ppi.p3().unit()),axis2.dot(ppi.p3().unit()));
	  if(phi<0.) phi+=2.*M_PI;
	  _h_piStar_phi->fill(phi);
	  // pions in lab frame
	  _h_pi_ctheta->fill(cos(q.polarAngle()));
	  _h_pi_phi   ->fill(q.azimuthalAngle());
	}
	// 2pi0
	else if (pi0.size()==2) {
	  FourMomentum q = pi0[0].momentum()+pi0[1].momentum();
	  _h_neutral_ex->fill(q.mass()/GeV);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_charged_ex);
      normalize(_h_charged_inc);
      normalize(_h_neutral_ex);
      normalize( _h_lept_ctheta);
      normalize(_h_lept_phi);
      normalize(_h_piStar_ctheta);
      normalize(_h_piStar_phi);
      normalize(_h_pi_ctheta);
      normalize(_h_pi_phi);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_charged_ex,_h_charged_inc,_h_neutral_ex;
    Histo1DPtr _h_lept_ctheta,_h_lept_phi,_h_piStar_ctheta,_h_piStar_phi,_h_pi_ctheta,_h_pi_phi;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1998_I467642);

}
