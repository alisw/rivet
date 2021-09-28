// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief Transverse Lambda polarization
  class BELLE_2019_I1687566 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2019_I1687566);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      declare(Thrust(fs),"Thrust");
      declare(UnstableParticles(),"UFS");
      declare(InitialQuarks(), "IQF");
      
      // Book histograms
      for(unsigned int ix=0;ix<4;++ix) {
      	book(_p_lam[ix]    ,1,ix+1,1);
      	book(_p_bar[ix]    ,1,ix+1,2);
      	book(_p_lam_pip[ix],2,ix+1,1);
      	book(_p_lam_pim[ix],2,ix+1,2);
      	book(_p_lam_Kp [ix],2,ix+1,3);
      	book(_p_lam_Km [ix],2,ix+1,4);
      	book(_p_bar_pip[ix],2,ix+1,5);
      	book(_p_bar_pim[ix],2,ix+1,6);
      	book(_p_bar_Kp [ix],2,ix+1,7);
      	book(_p_bar_Km [ix],2,ix+1,8);
      }
      book(_p_lam_inc   ,3,1,1);
      book(_p_lam_prompt,3,1,3);
      book(_p_lam_sigma ,3,1,5);
      book(_p_bar_inc   ,3,1,2);
      book(_p_bar_prompt,3,1,4);
      book(_p_bar_sigma ,3,1,6);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double alpha=0.642;
      static const double norm=3./alpha*100.;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;

      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      // thrust, to define an axis
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      const FinalState & fs = apply<FinalState>(event, "FS");
      
      if(thrust.thrust()<0.8)
	vetoEvent;
      for(const Particle & lambda : ufs.particles(Cuts::abspid==3122)) {
	int sign = lambda.pid()/3122;
	if(lambda.children().size()!=2) continue;
	// look at the decay products
	Particle proton,pion;
	if(lambda.children()[0].pid()==sign*2212 && 
	   lambda.children()[1].pid()==-sign*211) {
	  proton = lambda.children()[0];
	  pion   = lambda.children()[1];
	}
	else if(lambda.children()[1].pid()==sign*2212 && 
		lambda.children()[0].pid()==-sign*211) {
	  proton = lambda.children()[1];
	  pion   = lambda.children()[0];
	}
	else
	  continue;
	double xE = lambda.momentum().t()/meanBeamMom;
	if(xE<0.2 || xE>0.9) continue;
	Vector3 axis1 = lambda.momentum().p3().unit();
	// transverse polarization
	Vector3 axis2;
	if(lambda.momentum().p3().dot(thrust.thrustAxis())>=0.) {
	  axis2 = thrust.thrustAxis();
	}
	else {
	  axis2 =-thrust.thrustAxis();
	}
	Vector3 axis3 = axis2.cross(axis1).unit();	
	// boost to the lambda rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(lambda.momentum().betaVec());
	FourMomentum pproton = boost.transform(proton.momentum());
	int ibin=-1;
	if(xE<=0.3)
	  ibin=0;
	else if(xE<=0.4)
	  ibin=1;
	else if(xE<=0.5)
	  ibin=2;
	else if(xE<=0.9)
	  ibin=3;
	double cTheta = axis3.dot(pproton.p3().unit())*norm;
	double pT = sqrt(sqr(thrust.thrustMajorAxis().dot(lambda.momentum().p3()))+
			 sqr(thrust.thrustMinorAxis().dot(lambda.momentum().p3())));
	bool prompt=true;
	for(const Particle & parent : lambda.parents()) {
	  if(parent.abspid()==3212) {
	    prompt=false;
	    break;
	  }
	}
	// using thrust axis
	if(sign>0) {
	  _p_lam_inc->fill(xE,cTheta);
	  _p_lam[ibin]->fill(pT,cTheta);
	  if(flavour<=3) {
	    if(prompt) _p_lam_prompt->fill(xE,cTheta);
	    else       _p_lam_sigma ->fill(xE,cTheta);
	  }
	}
	else if(sign<0) {
	  _p_bar_inc->fill(xE,cTheta);
	  _p_bar[ibin]->fill(pT,cTheta);
	  if(flavour<=3) {
	    if(prompt) _p_bar_prompt->fill(xE,cTheta);
	    else       _p_bar_sigma ->fill(xE,cTheta);
	  }
	}
	for(const Particle & p2 : fs.particles(Cuts::abspid==211 || Cuts::abspid==321) ) {
	  if(axis2.dot(p2.p3())>0.) continue;
	  double xH = p2.momentum().t()/meanBeamMom;
	  Vector3 axisH = p2.momentum().p3();
	  Vector3 axis4 = axisH.cross(axis1).unit();
	  double cTheta2 = axis3.dot(pproton.p3().unit())*norm;
	  if(sign>0) {
	    if(p2.pid()==211)       _p_lam_pip[ibin]->fill(xH,cTheta2);
	    else if(p2.pid()==-211) _p_lam_pim[ibin]->fill(xH,cTheta2);
	    else if(p2.pid()== 321) _p_lam_Kp [ibin]->fill(xH,cTheta2);
	    else if(p2.pid()==-321) _p_lam_Km [ibin]->fill(xH,cTheta2);
	  }
	  else {
	    if(p2.pid()==211)       _p_bar_pip[ibin]->fill(xH,cTheta2);
	    else if(p2.pid()==-211) _p_bar_pim[ibin]->fill(xH,cTheta2);
	    else if(p2.pid()== 321) _p_bar_Kp [ibin]->fill(xH,cTheta2);
	    else if(p2.pid()==-321) _p_bar_Km [ibin]->fill(xH,cTheta2);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    //@}


    /// @name Histograms
    //@{
    Profile1DPtr _p_lam[4],_p_bar[4];
    Profile1DPtr _p_lam_pip[4],_p_lam_pim[4],_p_lam_Kp[4],_p_lam_Km[4];
    Profile1DPtr _p_bar_pip[4],_p_bar_pim[4],_p_bar_Kp[4],_p_bar_Km[4];
    Profile1DPtr _p_lam_inc,_p_lam_prompt,_p_lam_sigma;
    Profile1DPtr _p_bar_inc,_p_bar_prompt,_p_bar_sigma;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2019_I1687566);


}
