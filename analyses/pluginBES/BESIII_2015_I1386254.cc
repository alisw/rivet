// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> omega e+ nu_e
  class BESIII_2015_I1386254 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1386254);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==411);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      DD.addStable(PID::ETA);
      DD.addStable(PID::ETAPRIME);
      declare(DD, "DD");
      
      // Book histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, { 111,1}, {-11,1}, { 12,1}};
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	if(!DD.modeMatches(ix,5,mode)) continue;
       	const Particle & pip= DD.decayProducts()[ix].at( 211)[0];
       	const Particle & pim= DD.decayProducts()[ix].at(-211)[0];
       	const Particle & pi0= DD.decayProducts()[ix].at( 111)[0];
	const Particle & ep = DD.decayProducts()[ix].at(-11)[0];
	const Particle & nue= DD.decayProducts()[ix].at( 12)[0];
	FourMomentum pOmega = pip.momentum()+pim.momentum()+pi0.momentum(); 
        _h[0]->fill(pOmega.mass2());
        FourMomentum qq = DD.decaying()[ix].momentum()-pOmega;
        _h[1]->fill(qq.mass2());
	// boost momenta to D rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DD.decaying()[ix].momentum().betaVec());
       	FourMomentum pPP = boost.transform(pOmega);
      	Matrix3 ptoz(-pPP.p3().unit(), Vector3(0,0,1));
      	boost.preMult(ptoz);
      	// the momenta in frane to W along z
      	FourMomentum pD  = boost.transform(DD.decaying()[ix].momentum());
       	FourMomentum ppip = boost.transform(pip.momentum());
       	FourMomentum ppim = boost.transform(pim.momentum());
       	FourMomentum ppi0 = boost.transform(pi0.momentum());
      	FourMomentum pe  = boost.transform(ep .momentum());
      	FourMomentum pnu = boost.transform(nue.momentum());
       	pOmega = ppip+ppim+ppi0;
	qq = pD-pOmega;
	LorentzTransform boostOmega = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
	Vector3 n1 = (boostOmega.transform(ppip).p3().cross(boostOmega.transform(ppim).p3())).unit();
	_h[2]->fill(n1.dot(pOmega.p3().unit()));
       	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
       	Vector3 axisE = boostW.transform(pe).p3().unit();
       	_h[3]->fill(axisE.dot(qq.p3().unit()));
      // 	axisOmega.setZ(0.);
      // 	axisE.setZ(0.);
      // 	double chi = atan2(axisE.cross(axisOmega).dot(qq.p3().unit()), axisE.dot(axisOmega));
      // 	_h[imode+4]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1386254);

}
