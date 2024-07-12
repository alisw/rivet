// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D -> rho e nu
  class CLEO_2013_I1081165 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2013_I1081165);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==411 ||Cuts::pid==421);
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
      static const map<PdgId,unsigned int> & mode1 = { { 111,1}, {-211,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode2 = { { 211,1}, {-211,1}, {-11,1}, { 12,1}};
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	Particle pi2;
	if     (DD.decaying()[ix].pid()==421 && DD.modeMatches(ix,4,mode1)) {
	  pi2= DD.decayProducts()[ix].at(111)[0];
	}
      	else if(DD.decaying()[ix].pid()==411 && DD.modeMatches(ix,4,mode2)) {
	  pi2= DD.decayProducts()[ix].at(211)[0];
	}
	else continue;
	const Particle & pim= DD.decayProducts()[ix].at(-211)[0];
       	const Particle & ep = DD.decayProducts()[ix].at(-11)[0];
       	const Particle & nue= DD.decayProducts()[ix].at( 12)[0];
	FourMomentum pRho = pi2.momentum()+pim.momentum();
	if(abs(pRho.mass()-0.7753)>0.15) continue;
        FourMomentum qq = DD.decaying()[ix].momentum()-pRho;
        _h[0]->fill(qq.mass2());
      	// boost momenta to D rest frame
       	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DD.decaying()[ix].momentum().betaVec());
       	FourMomentum pPP = boost.transform(pRho);
      	Matrix3 ptoz(-pPP.p3().unit(), Vector3(0,0,1));
      	boost.preMult(ptoz);
       	// the momenta in frane to W along z
       	FourMomentum pD  = boost.transform(DD.decaying()[ix].momentum());
       	FourMomentum ppi2 = boost.transform(pi2.momentum());
       	FourMomentum ppim = boost.transform(pim.momentum());
      	FourMomentum pe  = boost.transform(ep .momentum());
      	FourMomentum pnu = boost.transform(nue.momentum());
       	pRho = ppi2+ppim;
       	qq = pD-pRho;
       	LorentzTransform boostRho = LorentzTransform::mkFrameTransformFromBeta(pRho.betaVec());
       	Vector3 axisRho = boostRho.transform(ppim).p3().unit();
	_h[1]->fill(axisRho.dot(pRho.p3().unit()));
      	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
	Vector3 axisE = boostW.transform(pe).p3().unit();
	_h[2]->fill(axisE.dot(qq.p3().unit()));
	axisRho.setZ(0.);
	axisE.setZ(0.);
	double chi = atan2(axisE.cross(axisRho).dot(qq.p3().unit()), axisE.dot(axisRho));
	_h[3]->fill(chi);
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


  RIVET_DECLARE_PLUGIN(CLEO_2013_I1081165);

}
