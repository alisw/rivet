// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- > D+ pi-pi-
  class BABAR_2009_I810694 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I810694);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 411);
      BP.addStable(-411);
      declare(BP, "BP");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix<2) book(_h_angle[ix],2,1,1+ix);
	book(_h_mass[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { -411,1},{ 211,2}};
      static const map<PdgId,unsigned int> & modeCC = { {  411,1},{-211,2}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	int sign=1;
      	if      (BP.modeMatches(ix,3,mode  )) sign= 1;
      	else if (BP.modeMatches(ix,3,modeCC)) sign=-1;
	else continue;
       	const Particle  & Dp  = BP.decayProducts()[ix].at(-sign*411)[0];
	const Particles & pim = BP.decayProducts()[ix].at( sign*211);
	_h_mass[2]->fill((pim[0].momentum()+pim[1].momentum()).mass2());
	// boost to B rest frame
	LorentzTransform boost =
	  LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. momentum().betaVec());
	FourMomentum pD    = boost.transform(Dp.momentum());
	FourMomentum ppi[2] = {boost.transform(pim[0].momentum()),boost.transform(pim[1].momentum())};
	double m2Dpi[2];
	for(unsigned int ix=0;ix<2;++ix) {
	  m2Dpi[ix] = (pim[ix].momentum()+Dp.momentum()).mass2();
	  if( (m2Dpi[ix]>4.5 && m2Dpi[ix]<5.5) ||
	      (m2Dpi[ix]>5.9 && m2Dpi[ix]<6.2) ) {
	    FourMomentum pDpi = pD+pim[ix];
	    LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDpi.betaVec());
	    Vector3 axis1 = boost2.transform(ppi[ix]).p3().unit();
	    Vector3 axis2 = (ix==0 ? ppi[1] : ppi[0]).p3().unit();
	    double cTheta = axis1.dot(axis2);
	    if(m2Dpi[ix]<5.5)
	      _h_angle[0]->fill(cTheta);
	    else
	      _h_angle[1]->fill(cTheta);
	  }
	}
	if(m2Dpi[0]>m2Dpi[1]) {
	  _h_mass[1]->fill(m2Dpi[0]);
	  _h_mass[0]->fill(m2Dpi[1]);
	}
	else {
	  _h_mass[0]->fill(m2Dpi[0]);
	  _h_mass[1]->fill(m2Dpi[1]);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix<2) normalize(_h_angle[ix],1.,false);
	normalize(_h_mass[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_angle[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I810694);

}
