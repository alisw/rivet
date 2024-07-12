// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> lambdabar p pi-
  class BABAR_2009_I819092 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I819092);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 3122);
      B0.addStable(-3122);
      declare(B0, "B0");
      book(_h_pol1,2,1,1);
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix<2) book(_h_mass[ix],1,1,1+ix);
	book(_h_pol2[ix],3,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1},{-3122,1}, {-211,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	if (!B0.modeMatches(ix,3,mode)) continue;
       	const Particle & pp     = B0.decayProducts()[ix].at( 2212)[0];
       	const Particle & LamBar = B0.decayProducts()[ix].at(-3122)[0];
	_h_mass[0]->fill( (pp.momentum()+LamBar.momentum()).mass());
	// boost to B rest frame
	LorentzTransform boost =
	  LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix]. momentum().betaVec());
	FourMomentum pLam    = boost.transform(LamBar.momentum());
	FourMomentum pProton = boost.transform(pp    .momentum());
	_h_mass[1]->fill(pLam.E());
	// Lambda decay products
	if(LamBar.children().size()!=2) continue;
	Particle pbar;
	if(LamBar.children()[0].pid()==-2212 &&
	   LamBar.children()[1].pid()== 211) {
	  pbar = LamBar.children()[0];
	}
	else if(LamBar.children()[1].pid()==-2212 &&
		LamBar.children()[0].pid()== 211) {
	  pbar = LamBar.children()[1];
	}
	else
	  continue;
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pLam.betaVec());
	Vector3 axisP = boost2.transform(boost.transform(pbar.momentum())).p3().unit();
	Vector3 axis1 = pLam.p3().unit();
	double cTheta = axisP.dot(axis1);
	_h_pol1   ->fill(pLam.E(),3.*cTheta);
	_h_pol2[0]->fill(pLam.E(),3.*cTheta);
	Vector3 axis2 = pLam.p3().cross(pProton.p3()).unit();
	cTheta = axisP.dot(axis2);
	_h_pol2[1]->fill(pLam.E(),3.*cTheta);
	Vector3 axis3 = axis1.cross(axis2);
	cTheta = axisP.dot(axis3);
	_h_pol2[2]->fill(pLam.E(),3.*cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double alpha = -0.732;
      for(unsigned int ix=0;ix<2;++ix) {
	if(ix<2) normalize(_h_mass[ix],1.,false);
	_h_pol2[ix]->scaleY(1./alpha);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    Profile1DPtr _h_pol1,_h_pol2[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I819092);

}
