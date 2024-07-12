// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K- pi+ e+ nu_e
  class BESIII_2016_I1411645 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1411645);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      DP.addStable(PID::K0S);
      DP.addStable(PID::ETA);
      DP.addStable(PID::ETAPRIME);
      declare(DP, "DP");
      
      // Book histograms
      for(unsigned int ix=0;ix<6;++ix)
	book(_h[ix],3,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { -321,1}, { 211,1}, {-11,1}, { 12,1}};
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      // loop over particles
      for(unsigned int ix=0;ix<DP.decaying().size();++ix) {
	if ( !DP.modeMatches(ix,4,mode) ) continue;
       	const Particle & Km = DP.decayProducts()[ix].at(-321)[0];
       	const Particle & pip= DP.decayProducts()[ix].at( 211)[0];
       	const Particle & ep = DP.decayProducts()[ix].at( -11)[0];
       	const Particle & nue= DP.decayProducts()[ix].at(  12)[0];
        FourMomentum pKstar = Km.momentum()+pip.momentum(); 
        _h[0]->fill(pKstar.mass());
        _h[1]->fill(pKstar.mass());
	FourMomentum qq = DP.decaying()[ix].momentum()-pKstar;
	_h[2]->fill(qq.mass2());
	// boost momenta to DP rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DP.decaying()[ix].momentum().betaVec());
	FourMomentum pKS = boost.transform(pKstar);
	Matrix3 ptoz(-pKS.p3().unit(), Vector3(0,0,1));
	boost.preMult(ptoz);
	// the momenta in frane to W along z
	FourMomentum pD  = boost.transform(DP.decaying()[ix].momentum());
	FourMomentum pK  = boost.transform(Km .momentum());
	FourMomentum ppi = boost.transform(pip.momentum());
	FourMomentum pe  = boost.transform(ep .momentum());
	FourMomentum pnu = boost.transform(nue.momentum());
	pKstar = pK+ppi;
	qq = pD-pKstar;
	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
      	Vector3 axisK = boostK.transform(pK).p3().unit();
      	_h[4]->fill(axisK.dot(pKstar.p3().unit()));
      	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
      	Vector3 axisE = boostW.transform(pe).p3().unit();
      	_h[3]->fill(axisE.dot(qq.p3().unit()));
      	axisK.setZ(0.);
      	axisE.setZ(0.);
      	double chi = atan2(axisE.cross(axisK).dot(qq.p3().unit()), axisE.dot(axisK));
      	_h[5]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1411645);

}
