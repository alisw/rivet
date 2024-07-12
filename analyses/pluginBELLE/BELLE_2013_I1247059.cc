// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> phi K* 
  class BELLE_2013_I1247059 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1247059);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(333);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 333,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 333,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
      	const Particle & Kp  = B0.decayProducts()[ix].at( 321*sign)[0];
      	const Particle & pim = B0.decayProducts()[ix].at(-211*sign)[0];
	const Particle & phi = B0.decayProducts()[ix].at( 333     )[0];
	if(phi.children().size()!=2 || phi.children()[0].pid()!=-phi.children()[1].pid() ||
	   phi.children()[0].abspid()!=321) continue;
	Particle Kp1 = phi.children()[0];
	Particle Km1 = phi.children()[1];
	if(Kp1.pid()<0) swap(Kp1,Km1);
	// B0 frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pKstar = boost1.transform(Kp.momentum()+pim.momentum());
	FourMomentum pPhi   = boost1.transform(phi.momentum());
	// stuff in K* frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	FourMomentum pKp = boost2.transform(boost1.transform(Kp.momentum()));
	Vector3 axis1 = pKstar.p3().unit();
	double cTheta1 = axis1.dot(pKp.p3().unit());
	if(cTheta1>0.75) continue;
	Vector3 trans1 = pKp.p3() - cTheta1*pKp.p3().mod()*axis1;
	// stuff in phi frame
	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
	FourMomentum pKp1 = boost3.transform(boost1.transform(Kp1.momentum()));
	Vector3 axis2 = pPhi.p3().unit();
	double cTheta2 = axis2.dot(pKp1.p3().unit());
	Vector3 trans2 = pKp1.p3() - cTheta2*pKp1.p3().mod()*axis2;
	// angle between planes
	double chi = atan2(trans1.cross(trans2).dot(axis1),trans1.dot(trans2));
	// fill histos
	_h[0]->fill(pKstar.mass());
	_h[1]->fill(cTheta1);
	_h[2]->fill(cTheta2);
	_h[3]->fill(chi);
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
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1247059);

}
