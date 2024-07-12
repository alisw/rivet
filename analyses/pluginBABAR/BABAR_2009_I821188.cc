// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 > KS0 pi+ pi-
  class BABAR_2009_I821188 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I821188);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::K0S);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 310,1}, { 211,1}, {-211,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	if (!B0.modeMatches(ix,3,mode)) continue;
      	int sign = B0.decaying()[ix].pid()>0 ? 1 : -1;
	// boost to B rest frame
	LorentzTransform boost =
	  LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix]. momentum().betaVec());
	// momenta
	FourMomentum pip  = boost.transform(B0.decayProducts()[ix].at( 211*sign)[0].momentum());
	FourMomentum pim  = boost.transform(B0.decayProducts()[ix].at(-211*sign)[0].momentum());
	FourMomentum K0   = boost.transform(B0.decayProducts()[ix].at( 310     )[0].momentum());
	// pi+pi- resonance
	FourMomentum ppipi = pim+pip;
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
	double cTheta = boost2.transform(pim).p3().unit().dot(K0.p3().unit());
	if(cTheta>0.) _h[0][0]->fill(ppipi.mass());
	else          _h[0][1]->fill(ppipi.mass());
	// K pi- resonance
	FourMomentum pKpim = K0+pim;
	boost2 = LorentzTransform::mkFrameTransformFromBeta(pKpim.betaVec());
	cTheta = boost2.transform(K0).p3().unit().dot(pip.p3().unit());
	if(cTheta>0.) _h[1][0]->fill(pKpim.mass());
	else          _h[1][1]->fill(pKpim.mass());
	// K pi+ resonance
	FourMomentum pKpip = K0+pip;
	boost2 = LorentzTransform::mkFrameTransformFromBeta(pKpip.betaVec());
	cTheta = boost2.transform(pip).p3().unit().dot(pim.p3().unit());
	if(cTheta>0.) _h[1][0]->fill(pKpip.mass());
	else          _h[1][1]->fill(pKpip.mass());
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


  RIVET_DECLARE_PLUGIN(BABAR_2009_I821188);

}
