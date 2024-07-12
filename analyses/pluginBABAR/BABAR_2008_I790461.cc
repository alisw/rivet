// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> K+K- e+ nu_e
  class BABAR_2008_I790461 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I790461);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      DS.addStable(PID::ETA);
      DS.addStable(PID::ETAPRIME);
      declare(DS, "DS");
      
      // Book histograms
      book(_h[0],1,1,1);
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix+1],2,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 321,1}, {-321,1}, {-11,1}, { 12,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	if(!DS.modeMatches(ix,4,mode)) continue;
	const Particle & Kp = DS.decayProducts()[ix].at( 321)[0];
	const Particle & Km = DS.decayProducts()[ix].at(-321)[0];
	const Particle & ep = DS.decayProducts()[ix].at( -11)[0];
	const Particle & nue= DS.decayProducts()[ix].at(  12)[0];
	FourMomentum pPhi = Kp.momentum()+Km.momentum(); 
	_h[0]->fill(pPhi.mass());
	FourMomentum qq = DS.decaying()[ix].momentum()-pPhi;
	_h[1]->fill(qq.mass2());
	// boost momenta to DS rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DS.decaying()[ix].momentum().betaVec());
	FourMomentum pPHI = boost.transform(pPhi);
	Matrix3 ptoz(-pPHI.p3().unit(), Vector3(0,0,1));
	boost.preMult(ptoz);
	// the momenta in frane to W along z
	FourMomentum pD  = boost.transform(DS.decaying()[ix].momentum());
	FourMomentum pKp = boost.transform(Kp .momentum());
	FourMomentum pKm = boost.transform(Km .momentum());
	FourMomentum pe  = boost.transform(ep .momentum());
	FourMomentum pnu = boost.transform(nue.momentum());
	pPhi = pKp+pKm;
	qq = pD-pPhi;
	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
	Vector3 axisK = boostK.transform(pKp).p3().unit();
	_h[3]->fill(axisK.dot(pPhi.p3().unit()));
	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(  qq.betaVec());
	Vector3 axisE = boostW.transform(pe).p3().unit();
	_h[2]->fill(axisE.dot(qq.p3().unit()));
	axisK.setZ(0.);
	axisE.setZ(0.);
	double chi = atan2(axisE.cross(axisK).dot(qq.p3().unit()), axisE.dot(axisK));
	_h[4]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I790461);

}
