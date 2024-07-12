// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma omega omega
  class BESIII_2013_I1203841 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1203841);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PHI);
      PSI.addStable(PID::OMEGA);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // histograms
      for(unsigned int ix=0;ix<9;++ix)
	book(_h[ix],1,1,1+ix);
    }

    // angle cuts due regions of BES calorimeter
    bool vetoPhoton(const double & cTheta) {
      return cTheta>0.92 || (cTheta>0.8 && cTheta<0.86);
    }
    
    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & pi0, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PI0) {
	  pi0.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,pi0,ncount);
      }
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
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 223,1}, { 333,1},{ 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particle & phi   = PSI.decayProducts()[0].at(333)[0];
      const Particle & omega = PSI.decayProducts()[0].at(223)[0];
      const Particle & gam   = PSI.decayProducts()[0].at( 22)[0];
      _h[0]->fill((omega.momentum()+phi.momentum()).mass());
      _h[1]->fill((omega.momentum()+gam.momentum()).mass());
      _h[2]->fill((phi  .momentum()+gam.momentum()).mass());
      double cTheta = axis.dot(gam.p3().unit());
      if(vetoPhoton(abs(cTheta))) vetoEvent;
      _h[3]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
      FourMomentum pGamma    = boost1.transform(gam.momentum());
      FourMomentum pOmegaPhi = boost1.transform(omega.momentum()+phi.momentum());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pOmegaPhi.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(phi.momentum())).p3().unit();
      _h[5]->fill(e1z.dot(axis2));
      double phiPhi = atan2(axis2.dot(e1y),axis2.dot(e1x));
      if(phiPhi<0.) phiPhi+=2.*M_PI;
      _h[7]->fill(phiPhi);
      // now for the phi decays
      if(phi.children().size()!=2|| phi.children()[0].pid()!=-phi.children()[1].pid() ||
	 phi.children()[0].abspid()!=321) vetoEvent;
      Particle Km = phi.children()[0];
      Particle Kp = phi.children()[1];
      if(Kp.pid()<0) swap(Km,Kp);
      FourMomentum pKp  = boost2.transform(boost1.transform(Kp.momentum()));
      FourMomentum pPhi = boost2.transform(boost1.transform(phi.momentum()));
      LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
      pKp = boost3.transform(pKp);
      double cK = axis2.dot(pKp.p3().unit());
      _h[6]->fill(cK);
      // omega decay
      unsigned int ncount=0;
      Particles pip,pim,pi0;
      findChildren(omega,pim,pip,pi0,ncount);
      if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) vetoEvent;
      // boost to omega/phi frame
      FourMomentum ppip = boost2.transform(boost1.transform(pip[0].momentum()));
      FourMomentum ppim = boost2.transform(boost1.transform(pim[0].momentum()));
      FourMomentum pOmega = boost2.transform(boost1.transform(omega.momentum()));
      LorentzTransform boost4 = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
      Vector3 axisZ = pOmega.p3().unit();
      ppip = boost4.transform(ppip);
      ppim = boost4.transform(ppim);
      Vector3 norm = ppip.p3().cross(ppim.p3()).unit();
      double cOmega = norm.dot(axisZ);
      _h[4]->fill(cOmega);
      // angle between planes
      Vector3 Trans1 = pKp.p3() - cK*pKp.p3().mod()*axis2;
      Vector3 Trans2 = norm     - cOmega*axisZ;
      double chi = atan(Trans1.cross(Trans2).dot(axis2)/Trans1.dot(Trans2));
      _h[8]->fill(abs(chi)/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1203841);

}
