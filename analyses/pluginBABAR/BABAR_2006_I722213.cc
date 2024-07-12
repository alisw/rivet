// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> Lambda Kbar0 K+
  class BABAR_2006_I722213 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I722213);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable(PID::ETA);
      LAMBDAC.addStable(PID::LAMBDA);
      LAMBDAC.addStable(-PID::LAMBDA);
      declare(LAMBDAC, "LAMBDAC");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_mass[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_ctheta[ix],2,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::LAMBDA,1}, { 321,1}, { 310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-PID::LAMBDA,1}, {-321,1}, { 310,1}};
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for(unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
	int sign = 1;
	if (LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & lam = LAMBDAC.decayProducts()[ix].at( sign*PID::LAMBDA)[0];
	const Particle & K0  = LAMBDAC.decayProducts()[ix].at( 310)[0];
	const Particle & Kp  = LAMBDAC.decayProducts()[ix].at( sign*321)[0];
	double mLamK = (lam.momentum()+K0.momentum()).mass();
	_h_mass[0]->fill(mLamK);
	_h_mass[1]->fill((Kp .momentum()+K0.momentum()).mass()-Kp.mass()-K0.mass());
	_h_mass[2]->fill((lam.momentum()+Kp.momentum()).mass()-Kp.mass()-lam.mass());
	// take Xi(1690) to be any resonance in mass region
	if(mLamK<1.6725 || mLamK>1.6975) continue;
	if(LAMBDAC.decaying()[ix].children().size()!=2) continue;
	if(LAMBDAC.decaying()[ix].children()[0].abspid()!=321 &&
	   LAMBDAC.decaying()[ix].children()[1].abspid()!=321) continue;
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(LAMBDAC.decaying()[ix].momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(LAMBDAC.decaying()[ix].momentum());
	FourMomentum pbaryon2 = boost1.transform(lam.momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_ctheta[0]->fill(cTheta);
	_h_ctheta[1]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_mass  [ix]);
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_ctheta[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass  [3];
    Histo1DPtr _h_ctheta[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I722213);

}
