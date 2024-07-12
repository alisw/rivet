// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> phi phi K
  class BABAR_2011_I901433 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2011_I901433);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PHI);
      BB.addStable(PID::K0S);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<6;++ix)
	book(_h[ix],1,1,1+ix);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      static const map<PdgId,unsigned int> & mode1   = { { 333,2},{ 321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 333,2},{-321,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 333,2},{ 310,1}};
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	if(BB.modeMatches(ix,3,mode1) || BB.modeMatches(ix,3,mode1CC) ||
	   BB.modeMatches(ix,3,mode2)) {
	  // phi mesons
	  const Particles & phi = BB.decayProducts()[ix].at(333);
	  // bost to B rest frane
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	  FourMomentum pPhiPhi= boost1.transform(phi[0].momentum()+phi[1].momentum());
	  double mPhiPhi = pPhiPhi.mass();
	  int iloc=-1;
	  if(mPhiPhi>2.94 && mPhiPhi<3.02)
	    iloc=0;
	  else if(mPhiPhi<2.85)
	    iloc=3;
	  else continue;
	  // cos theta phi phi
	  Vector3 axis1 = pPhiPhi.p3().unit();
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pPhiPhi.betaVec());
	  Vector3 axis2 = boost2.transform(boost1.transform(phi[0].momentum())).p3().unit();
	  _h[iloc+2]->fill(abs(axis1.dot(axis2)));
	  // now for the phi decays
	  Vector3 Trans[2];
	  bool foundPhi=true;
	  for(unsigned int ix=0;ix<2;++ix) {
	    if(phi[ix].children().size()!=2|| phi[ix].children()[0].pid()!=-phi[ix].children()[1].pid() ||
	       phi[ix].children()[0].abspid()!=321) {
	      foundPhi = false;
	      break;
	    }
	    Particle Km = phi[ix].children()[0];
	    Particle Kp = phi[ix].children()[1];
	    if(Kp.pid()<0) swap(Km,Kp);
	    FourMomentum pKp  = boost2.transform(boost1.transform(Kp.momentum()));
	    FourMomentum pPhi = boost2.transform(boost1.transform(phi[ix].momentum()));
	    LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
	    pKp = boost3.transform(pKp);
	    double cK = axis2.dot(pKp.p3().unit());
	    _h[iloc+1]->fill(cK);
	    Trans[ix] = pKp.p3() - cK*pKp.p3().mod()*axis2;
	  }
	  if(!foundPhi) continue;
	  double chi = atan2(Trans[0].cross(Trans[1]).dot(axis2),Trans[0].dot(Trans[1]));
	  _h[iloc]->fill(abs(chi));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2011_I901433);

}
