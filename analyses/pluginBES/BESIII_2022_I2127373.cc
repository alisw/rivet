// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> Lambda0 e+ nu_e
  class BESIII_2022_I2127373 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2127373);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable(PID::ETA);
      LAMBDAC.addStable(PID::ETAPRIME);
      declare(LAMBDAC, "LAMBDAC");
      
      // Book histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 2212,1}, {-211,1}, {-11,1}, { 12,1}};
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for(unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
      	if ( !LAMBDAC.modeMatches(ix,4,mode) ) continue;
	const Particle & pp = LAMBDAC.decayProducts()[ix].at(2212)[0];
	const Particle & pim= LAMBDAC.decayProducts()[ix].at(-211)[0];
	const Particle & ep = LAMBDAC.decayProducts()[ix].at( -11)[0];
	const Particle & nue= LAMBDAC.decayProducts()[ix].at(  12)[0];
	if(LAMBDAC.decaying()[ix].children(Cuts::pid==PID::LAMBDA).empty()) continue;
	FourMomentum pLambda = pp.momentum()+pim.momentum(); 
	FourMomentum qq = LAMBDAC.decaying()[ix].momentum()-pLambda;
	_h[0]->fill(qq.mass2());
	// boost momenta to LAMBDAC rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(LAMBDAC.decaying()[ix].momentum().betaVec());
	FourMomentum pLam = boost.transform(pLambda);
	Matrix3 ptoz(-pLam.p3().unit(), Vector3(0,0,1));
	boost.preMult(ptoz);
	// the momenta in frane to W along z
	FourMomentum pD  = boost.transform(LAMBDAC.decaying()[ix].momentum());
	FourMomentum pP  = boost.transform(pp .momentum());
	FourMomentum ppi = boost.transform(pim.momentum());
	FourMomentum pe  = boost.transform(ep .momentum());
	FourMomentum pnu = boost.transform(nue.momentum());
	pLambda = pP+ppi;
	qq = pD-pLambda;
	LorentzTransform boostL = LorentzTransform::mkFrameTransformFromBeta(pLambda.betaVec());
	Vector3 axisP = boostL.transform(pP).p3().unit();
	_h[1]->fill(axisP.dot(pLambda.p3().unit()));
	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
	Vector3 axisE = boostW.transform(pe).p3().unit();
	_h[2]->fill(-axisE.dot(qq.p3().unit()));
	axisP.setZ(0.);
	axisE.setZ(0.);
	double chi = atan2(axisE.cross(axisP).dot(qq.p3().unit()), axisE.dot(axisP));
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


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2127373);

}
