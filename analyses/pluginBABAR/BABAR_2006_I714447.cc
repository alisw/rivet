// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D*(sJ) production
  class BABAR_2006_I714447 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2006_I714447);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      // dists
      book(_s_2317  ,2,1,1);
      book(_s_2460_1,3,1,1);
      book(_s_2460_2,4,1,1);
      book(_hel     ,5,1,1);
      // rates
      book(_r_2317  ,1,1,1);
      book(_r_2460_1,1,1,2);
      book(_r_2460_2,1,1,3);
      book(_r_2460_3,1,1,4);
      book(_r_2460_4,1,1,5);
      book(_r_2536  ,1,1,6);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & Ds: apply<UnstableParticles>(event,"UFS").particles(Cuts::pid==10433 or
									       Cuts::pid==20433 or
									       Cuts::pid==10431)) {
	// cut of 3.2 GeV of the momentum
	double mom=Ds.momentum().p3().mod();
	if(mom<=3.2) continue;
	if(Ds.abspid()==10431) {
	  _s_2317->fill(mom);
	}
	else if(Ds.abspid()==20433) {
	  _s_2460_1->fill(mom);
	  _s_2460_2->fill(mom);
	}
	int sign = Ds.pid()/Ds.abspid();
	Particle DsStar;
	bool doHelicity=false;
	// identify the decay mode
	if(Ds.children().size()==2) {
	  if(Ds.children()[0].pid()==sign*431 &&
	     Ds.children()[1].pid()==111) {
	    if(Ds.abspid()==10431) _r_2317->fill(10.58);
	  }
	  else if(Ds.children()[1].pid()==sign*431 &&
		  Ds.children()[0].pid()==111) {
	    if(Ds.abspid()==10431) _r_2317->fill(10.58);
	  }
	  else if(Ds.children()[0].pid()==sign*433 &&
		  Ds.children()[1].pid()==111) {
	    if(Ds.abspid()==20433) {
	      _r_2460_3->fill(10.58);
	      DsStar=Ds.children()[0];
	      doHelicity=true;
	    }
	  }
	  else if(Ds.children()[1].pid()==sign*433 &&
		  Ds.children()[0].pid()==111) {
	    if(Ds.abspid()==20433) {
	      _r_2460_3->fill(10.58);
	      DsStar=Ds.children()[1];
	      doHelicity=true;
	    }
	  }
	  else if(Ds.children()[0].pid()==sign*431 &&
		  Ds.children()[1].pid()==22) {
	    if(Ds.abspid()==20433) _r_2460_1->fill(10.58);
	  }
	  else if(Ds.children()[1].pid()==sign*431 &&
		  Ds.children()[0].pid()==22) {
	    if(Ds.abspid()==20433) _r_2460_1->fill(10.58);
	  }
	  if(doHelicity && DsStar.children().size()==2) {
	    Particle gamma;
	    doHelicity=false;
	    for(const Particle & p : DsStar.children()) {
	      if(p.pid()==22) {
		doHelicity=true;
		gamma=p;
	      }
	    }
	    if(doHelicity) {
	      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(Ds.momentum().betaVec());
	      FourMomentum pDsStar = boost.transform(DsStar.momentum());
	      FourMomentum pGamma  = boost.transform(gamma .momentum());
	      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDsStar.betaVec());
	      FourMomentum pGamma2  = boost2.transform(pGamma);
	      double cTheta = pGamma2.p3().unit().dot(pDsStar.p3().unit());
	      _hel->fill(cTheta);
	    }
	  }
	}
	Particles D, pip, pi0, gamma;
	unsigned int nstable(0);
	findDecayProducts(Ds,D,pip,pi0,gamma,nstable);
	if(nstable==3) {
	  if(D.size()==1&&pip.size()==2) {
	    if(Ds.abspid()==20433) {
	      _r_2460_4->fill(10.58);
	    }
	    else if(Ds.abspid()==10433) {
	      _r_2536  ->fill(10.58);
	    }
	  }
	  else if(D.size()==1&&pi0.size()==1&&gamma.size()==1) {
	    if(Ds.abspid()==20433)
	      _r_2460_2->fill(10.58);
	  }							  
	}
      }
    }

    void findDecayProducts(Particle parent, Particles & Ds, Particles & pip, Particles & pi0,
			   Particles & gamma, unsigned int & nstable) {
      for(const Particle & p: parent.children()) {
	if(p.pid()==111) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if(p.pid()==22) {
	  gamma.push_back(p);
	  ++nstable;
	}
	else if(p.abspid()==211) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if(p.abspid()==431) {
	  Ds.push_back(p);
	  ++nstable;
	}
	else if(p.children().empty()) {
	  ++nstable;
	}
	else
	  findDecayProducts(p,Ds,pip,pi0,gamma,nstable);
      }
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_s_2317  ,1.,false);
      normalize(_s_2460_1,1.,false);
      normalize(_s_2460_2,1.,false);
      normalize(_hel     ,1.,false);
      scale(_r_2317  ,crossSection()/femtobarn/sumOfWeights());
      scale(_r_2460_1,crossSection()/femtobarn/sumOfWeights());
      scale(_r_2460_2,crossSection()/femtobarn/sumOfWeights());
      scale(_r_2460_3,crossSection()/femtobarn/sumOfWeights());
      scale(_r_2460_4,crossSection()/femtobarn/sumOfWeights());
      scale(_r_2536  ,crossSection()/femtobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _s_2317,_s_2460_1,_s_2460_2,_hel;
    Histo1DPtr _r_2317,_r_2460_1,_r_2460_2,_r_2460_3,_r_2460_4,_r_2536;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2006_I714447);

}
