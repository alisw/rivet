// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> Kbar0 pi- e+ nu_e
  class BESIII_2018_I1705754 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1705754);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      
      // Book histograms
      for(unsigned int ix=0;ix<5;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 310,1}, {-211,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode2 = { { 130,1}, {-211,1}, {-11,1}, { 12,1}};
      static const map<PdgId,unsigned int> & mode3 = { {-311,1}, {-211,1}, {-11,1}, { 12,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	Particle K0;
       	if     (D0.modeMatches(ix,4,mode1)) K0=D0.decayProducts()[ix].at( 310)[0];
	else if(D0.modeMatches(ix,4,mode2)) K0=D0.decayProducts()[ix].at( 130)[0];
	else if(D0.modeMatches(ix,4,mode3)) K0=D0.decayProducts()[ix].at(-311)[0];
	else continue;
	const Particle & pim= D0.decayProducts()[ix].at(-211)[0];
	const Particle & ep = D0.decayProducts()[ix].at(-11)[0];
	const Particle & nue= D0.decayProducts()[ix].at( 12)[0];
        FourMomentum pKstar = K0.momentum()+pim.momentum(); 
        _h[0]->fill(pKstar.mass());
        FourMomentum qq = D0.decaying()[ix].momentum()-pKstar;
        _h[1]->fill(qq.mass2());
	// boost momenta to D0 rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(D0.decaying()[ix].momentum().betaVec());
	FourMomentum pKS = boost.transform(pKstar);
	Matrix3 ptoz(-pKS.p3().unit(), Vector3(0,0,1));
	boost.preMult(ptoz);
	// the momenta in frane to W along z
	FourMomentum pD  = boost.transform(D0.decaying()[ix].momentum());
	FourMomentum pK  = boost.transform(K0 .momentum());
	FourMomentum ppi = boost.transform(pim.momentum());
	FourMomentum pe  = boost.transform(ep .momentum());
	FourMomentum pnu = boost.transform(nue.momentum());
	pKstar = pK+ppi;
	qq = pD-pKstar;
	LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	Vector3 axisK = boostK.transform(pK).p3().unit();
	_h[3]->fill(axisK.dot(pKstar.p3().unit()));
	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
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


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1705754);

}
