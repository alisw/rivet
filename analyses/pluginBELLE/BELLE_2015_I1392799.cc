// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> p Lambda D
  class BELLE_2015_I1392799 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1392799);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 411);
      B0.addStable(-411);
      B0.addStable( 413);
      B0.addStable(-413);
      B0.addStable(PID::LAMBDA);
      B0.addStable(-PID::LAMBDA);
      declare(B0, "B0");
      //histograms
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
      book(_nB,"/TMP/nb");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 2212,1},{-3122,1}, {-411,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-2212,1},{ 3122,1}, { 411,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 2212,1},{-3122,1}, {-413,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-2212,1},{ 3122,1}, { 413,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	_nB->fill();
      	int sign = 1, imode = 0;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode1)) {
	  imode=0;
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode1CC)) {
	  imode=0;
      	  sign=-1;
      	}
	else if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode2)) {
	  imode=1;
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode2CC)) {
	  imode=1;
      	  sign=-1;
      	}
      	else
      	  continue;
	int iD = imode==0 ? 411 : 413;
	const Particle & pp  = B0.decayProducts()[ix].at( sign*2212)[0];
	const Particle & lam = B0.decayProducts()[ix].at(-sign*3122)[0];
	const Particle & Dm  = B0.decayProducts()[ix].at(-sign*iD  )[0];
	_h[0][imode]->fill((lam.momentum()+pp.momentum()).mass());
	// first boost to B rest frame
	LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum ppp  = boostB.transform(pp .momentum());
	FourMomentum plam = boostB.transform(lam.momentum());
	FourMomentum pDm  = boostB.transform(Dm .momentum());
	LorentzTransform boostpLam = LorentzTransform::mkFrameTransformFromBeta((ppp+plam).betaVec());
	ppp = boostpLam.transform(ppp);
	pDm  = boostpLam.transform(pDm);
	_h[1][imode]->fill(ppp.p3().unit().dot(pDm.p3().unit()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	double fact = ix==0 ? 1e3 : 1e6;
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h[ix][iy], fact/ *_nB);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    CounterPtr _nB;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1392799);

}
