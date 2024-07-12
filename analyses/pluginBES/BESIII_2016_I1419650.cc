// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma phi phi
  class BESIII_2016_I1419650 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1419650);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PHI);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // book histograms
      for(unsigned int ix=0;ix<5;++ix)
	book(_h[ix],1,1,1+ix);
    }

    // angle cuts due regions of BES calorimeter
    bool vetoPhoton(const double & cTheta) {
      return cTheta>0.92 || (cTheta>0.8 && cTheta<0.86);
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
      static const map<PdgId,unsigned int> & mode = { { 333,2},{ 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particles & phi = PSI.decayProducts()[0].at(333);
      const Particle  & gam = PSI.decayProducts()[0].at( 22)[0];
      _h[0]->fill((phi[0].momentum()+phi[1].momentum()).mass());
      double cTheta = axis.dot(gam.p3().unit());
      if(vetoPhoton(abs(cTheta))) vetoEvent;
      _h[1]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
      FourMomentum pGamma = boost1.transform(gam.momentum());
      FourMomentum pPhiPhi= boost1.transform(phi[0].momentum()+phi[1].momentum());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pPhiPhi.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(phi[0].momentum())).p3().unit();
      _h[2]->fill(e1z.dot(axis2));
      // now for the phi decays
      Particle Km[2],Kp[2];
      FourMomentum pKp[2],pPhi[2];
      for(unsigned int ix=0;ix<2;++ix) {
	if(phi[ix].children().size()!=2|| phi[ix].children()[0].pid()!=-phi[ix].children()[1].pid() ||
	   phi[ix].children()[0].abspid()!=321) vetoEvent;
	Km[ix] = phi[ix].children()[0];
	Kp[ix] = phi[ix].children()[1];
	if(Kp[ix].pid()<0) swap(Km[ix],Kp[ix]);
	pKp[ix]  = boost2.transform(boost1.transform(Kp[ix].momentum()));
	pPhi[ix] = boost2.transform(boost1.transform(phi[ix].momentum()));
	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi[ix].betaVec());
	pKp[ix] = boost3.transform(pKp[ix]);
      }
      double cK = axis2.dot(pKp[0].p3().unit());
      _h[3]->fill(cK);
      Vector3 Trans1 = pKp[0].p3() - cK*pKp[0].p3().mod()*axis2;
      Vector3 Trans2 = pKp[1].p3() - axis2.dot(pKp[1].p3())*axis2;
      double chi = atan(Trans1.cross(Trans2).dot(axis2)/Trans1.dot(Trans2));
      _h[4]->fill(abs(chi)/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1419650);

}
