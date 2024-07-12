// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief q^2 in D_s+ -> K0 e+ nu_e and ->K*0 e+ nu_e distributions
  class BESIII_2019_I1702549 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1702549);


    /// @name Analysis methods
    ///@{

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
      book(_h_q2, 1, 1, 1);
      book(_nD,"/TMP/nD");
      for(unsigned int ix=0;ix<5;++ix)
	book(_h_Kstar[ix],1,1,3+ix);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 310,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode2 = { { 130,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode3 = { { 311,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode4 = { { 321,1}, {-211,1}, {-11,1}, { 12,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	_nD->fill();
	if(DS.modeMatches(ix,3,mode1)) {
	  _h_q2->fill((DS.decaying()[ix].momentum()-DS.decayProducts()[ix].at(310)[0].momentum()).mass2());
	}
	else if(DS.modeMatches(ix,3,mode2)) {
	  _h_q2->fill((DS.decaying()[ix].momentum()-DS.decayProducts()[ix].at(130)[0].momentum()).mass2());
	}
	else if(DS.modeMatches(ix,3,mode3)) {
	  _h_q2->fill((DS.decaying()[ix].momentum()-DS.decayProducts()[ix].at(311)[0].momentum()).mass2());
	}
	else if(DS.modeMatches(ix,4,mode4)) {
	  const Particle & Kp = DS.decayProducts()[ix].at( 321)[0];
	  const Particle & pim= DS.decayProducts()[ix].at(-211)[0];
	  const Particle & ep = DS.decayProducts()[ix].at(-11)[0];
	  const Particle & nue= DS.decayProducts()[ix].at( 12)[0];
	  FourMomentum pKstar = Kp.momentum()+pim.momentum(); 
	  _h_Kstar[0]->fill(pKstar.mass());
	  FourMomentum qq = DS.decaying()[ix].momentum()-pKstar;
	  _h_Kstar[1]->fill(qq.mass2());
	  // boost momenta to Ds rest frame
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DS.decaying()[ix].momentum().betaVec());
	  FourMomentum pKS = boost.transform(pKstar);
	  Matrix3 ptoz(-pKS.p3().unit(), Vector3(0,0,1));
	  boost.preMult(ptoz);
	  // the momenta in frane to W along z
	  FourMomentum pD  = boost.transform(DS.decaying()[ix].momentum());
	  FourMomentum pK  = boost.transform(Kp .momentum());
	  FourMomentum ppi = boost.transform(pim.momentum());
	  FourMomentum pe  = boost.transform(ep .momentum());
	  FourMomentum pnu = boost.transform(nue.momentum());
	  pKstar = pK+ppi;
	  qq = pD-pKstar;
	  LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	  Vector3 axisK = boostK.transform(pK).p3().unit();
	  _h_Kstar[3]->fill(axisK.dot(pKstar.p3().unit()));
	  LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
	  Vector3 axisE = boostW.transform(pe).p3().unit();
	  _h_Kstar[2]->fill(axisE.dot(qq.p3().unit()));
	  axisK.setZ(0.);
	  axisE.setZ(0.);
	  double chi = atan2(axisE.cross(axisK).dot(qq.p3().unit()), axisE.dot(axisK));
          _h_Kstar[4]->fill(chi);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       // normalise to width in inverse ns
       scale(_h_q2, 1./0.504e-3/ *_nD);
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h_Kstar[ix]);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_q2;
    Histo1DPtr _h_Kstar[5];
    CounterPtr _nD;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1702549);

}
