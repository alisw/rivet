// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K- pi+ e+ nu_e
  class BABAR_2010_I879997 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I879997);


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
      for(unsigned int ix=0;ix<5;++ix)
	book(_h[ix],1,1,1+ix);
      double bins[5] = {0.,0.8,0.9,1.,1.6};
      for(unsigned int ix=0;ix<4;++ix) {
	for(unsigned int iy=0;iy<4;++iy) {
	  Histo1DPtr tmp;
	  book(tmp,2+iy,1,1+ix);
	  _b[ix].add(bins[iy],bins[iy+1],tmp);
	}
      }
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
	double mKpi = pKstar.mass();
	_h[4]->fill(mKpi);
	FourMomentum qq = DP.decaying()[ix].momentum()-pKstar;
	double q2 = qq.mass2();
	_h[0]->fill(q2);
	_b[0].fill(mKpi,q2);
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
	double cosK = axisK.dot(pKstar.p3().unit());
      	_h[2]->fill(cosK);
      	_b[2].fill(mKpi,cosK);
      	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(    qq.betaVec());
      	Vector3 axisE = boostW.transform(pe).p3().unit();
	double cosE = axisE.dot(qq.p3().unit());
      	_h[3]->fill(cosE);
	_b[3].fill(mKpi,cosE);
      	axisK.setZ(0.);
      	axisE.setZ(0.);
      	double chi = atan2(axisE.cross(axisK).dot(qq.p3().unit()), axisE.dot(axisK));
      	_h[1]->fill(chi);
	_b[1].fill(mKpi,chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix)
	normalize(_h[ix]);
      for(unsigned int ix=0;ix<4;++ix)
	for (Histo1DPtr hist : _b[ix].histos())
	  normalize(hist);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    BinnedHistogram _b[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I879997);

}
