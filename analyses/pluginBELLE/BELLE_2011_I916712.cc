// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> K X(3872)
  class BELLE_2011_I916712 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2011_I916712);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(UnstableParticles(Cuts::abspid==511 ||
				Cuts::abspid==521), "UFS");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1+ix,1,1);
    }
    
    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & Jpsi, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::JPSI) {
	  Jpsi.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,Jpsi,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles()) {
	if(p.children().empty()) continue;
	if(p.children().size()==1) continue;
	if(p.children().size()!=2) continue;
	Particle K,X;
	if(p.children()[0].pid()==_pid) {
	  X = p.children()[0];
	  K = p.children()[1];
	}
	else if(p.children()[1].pid()==_pid) {
	  X = p.children()[1];
	  K = p.children()[0];
	}
	else continue;
	if(K.abspid()!=311 && K.abspid()!=321 &&
	   K.abspid()!=310 && K.abspid()!=130) continue;
	// X(3872) decay
	unsigned int ncount=0;
	Particles pip,pim,Jpsi;
	findChildren(X,pim,pip,Jpsi,ncount);
	if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && Jpsi.size()==1)) continue;
	_h[3]->fill((pip[0].momentum()+pim[0].momentum()).mass());
	LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axisX = -boostB.transform(K.momentum()).p3().unit();
	FourMomentum pX   = boostB.transform(X      .momentum());
	LorentzTransform boostX = LorentzTransform::mkFrameTransformFromBeta(pX.betaVec());
	FourMomentum pPsi = boostX.transform(boostB.transform(Jpsi[0].momentum()));
	double cTheta = axisX.dot(pPsi.p3().unit());
	_h[0]->fill(cTheta);
	// finally the leptons from J/psi decay
	if(Jpsi[0].children().size()!=2) vetoEvent;
	if(Jpsi[0].children()[0].pid()!=-Jpsi[0].children()[1].pid()) vetoEvent;
	if(Jpsi[0].children()[0].abspid()!=PID::EMINUS &&
	   Jpsi[0].children()[0].abspid()!=PID::MUON) vetoEvent;
	Particle lm = Jpsi[0].children()[0];
	Particle lp = Jpsi[0].children()[1];
	Vector3 axispi = boostX.transform(boostB(pip[0].momentum())).p3().unit();
	Vector3 axisZ = axispi.cross(axisX).unit();
	Vector3 axisL = boostX.transform(boostB(lp.momentum())).p3().unit();
	_h[2]->fill(abs(axisZ.dot(axisL)));
	_h[1]->fill(abs(axisX.dot(axisL)));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2011_I916712);

}
