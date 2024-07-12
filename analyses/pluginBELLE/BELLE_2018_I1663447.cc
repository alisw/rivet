// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> K+ eta gamma
  class BELLE_2018_I1663447 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1663447);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::ETA);
      declare(BP, "BP");
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 221,1},{ 321,1}, {22,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 221,1},{-321,1}, {22,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & Kp    = BP.decayProducts()[ix].at( sign*321)[0];
	const Particle & eta   = BP.decayProducts()[ix].at(      221)[0];
	const Particle & gamma = BP.decayProducts()[ix].at(       22)[0];
	FourMomentum pKeta = Kp.momentum()+eta.momentum();
	double mass = pKeta.mass();
	LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].momentum().betaVec());
	pKeta = boostB.transform(pKeta);
	LorentzTransform boostKeta = LorentzTransform::mkFrameTransformFromBeta(pKeta.betaVec());
	FourMomentum pK     = boostKeta.transform(boostB.transform(Kp   .momentum()));
	FourMomentum pGamma = boostB.transform(gamma.momentum());
	double cTheta = pK.p3().unit().dot(pGamma.p3().unit());
	for(unsigned int iy=0;iy<2;++iy) {
	  _h[0][iy]->fill(cTheta);
	  _h[1][iy]->fill(mass );
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1663447);

}
