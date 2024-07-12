// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> pi+ pi+ pi-
  class BABAR_2009_I813140 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I813140);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      declare(BP, "BP");
      // histograms
      book(_h_all,1,1,1);
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_charge[ix][iy],2,1+ix,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { -211,1},{ 211,2}};
      static const map<PdgId,unsigned int> & modeCC = { {  211,1},{-211,2}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	int sign=1;
	if      (BP.modeMatches(ix,3,mode  )) sign= 1;
	else if (BP.modeMatches(ix,3,modeCC)) sign=-1;
	else continue;
	// particles
	const Particle  & pim = BP.decayProducts()[ix].at(-sign*211)[0];
       	const Particles & pip = BP.decayProducts()[ix].at( sign*211);
     	// boost to B rest frame
     	LorentzTransform boost =
     	  LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. momentum().betaVec());
	FourMomentum pPim = boost.transform(pim.momentum());
	FourMomentum pPip[2];
	for(unsigned int ix=0;ix<2;++ix) pPip[ix] = boost.transform(pip[ix].momentum());
	// loop over pi+
	for(unsigned int ix=0;ix<2;++ix) {
	  FourMomentum ppipi = pip[ix].momentum()+pim.momentum(); 
	  double mpipi = ppipi.mass();
	  _h_all->fill(mpipi);
	  _h_charge[0][(1-sign)/2]->fill(mpipi);
	  unsigned int ibatch = ix==0 ? 1 : 0;
	  Vector3 axis = pPip[ibatch].p3().unit();
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
	  Vector3 axis2 = boost1.transform(pPim).p3().unit();
	  double cTheta = axis.dot(axis2);
	  if(cTheta>0) {
	    _h_charge[1][(1-sign)/2]->fill(mpipi);
	  }
	  else {
	    _h_charge[2][(1-sign)/2]->fill(mpipi);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_all,1.,false);
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h_charge[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_all,_h_charge[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I813140);

}
