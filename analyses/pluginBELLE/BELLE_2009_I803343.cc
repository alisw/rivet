// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> Lambda Lambdabar K(*)0
  class BELLE_2009_I803343 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I803343);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==511 ||
						Cuts::pid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 3122);
      BB.addStable(-3122);
      BB.addStable( 310);
      BB.addStable( 313);
      BB.addStable(-313);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_mass[ix],1,1,1+ix);
	book(_h_angle[ix],2+ix,1,1);
      }
      book(_c[0],"TMP/nB0");
      book(_c[1],"TMP/nBP");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 3122,1},{-3122,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 3122,1},{-3122,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 3122,1},{-3122,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 3122,1},{-3122,1}, { 313,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 3122,1},{-3122,1}, {-313,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	if(BB.decaying()[ix].abspid()==511) _c[0]->fill();
	else                                _c[1]->fill();
	unsigned int imode=0;
	if (BB.modeMatches(ix,3,mode1))
	  imode=0;
	else if((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) ||
		(BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)))
	  imode=1;
	else if((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) ||
		(BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC)))
	  imode=2;
	else
	  continue;
	int sign = BB.decaying()[ix].pid()>0 ? 1 : -1;
       	const Particle & Lam    = BB.decayProducts()[ix].at( sign*3122)[0];
       	const Particle & LamBar = BB.decayProducts()[ix].at(-sign*3122)[0];
	FourMomentum pLL = Lam.momentum()+LamBar.momentum();
	double mass = pLL.mass();
	_h_mass[imode]->fill(mass);
	// rest just in threshold region
	if(mass>2.85) continue;
	// boost to B rest frame
	LorentzTransform boost =
	  LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix]. momentum().betaVec());
	// B+ K+
	if(imode==1) {
	  pLL = boost.transform(pLL);
	  LorentzTransform boost2 =
	    LorentzTransform::mkFrameTransformFromBeta(pLL.betaVec());
	  
	  FourMomentum pLam    = boost2.transform(boost.transform(Lam.momentum()));
	  FourMomentum pLamB   = boost2.transform(boost.transform(LamBar.momentum()));
	  const Particle & Kp  = BB.decayProducts()[ix].at( sign*321)[0];
	  FourMomentum pK      = boost2.transform(boost.transform(Kp.momentum()));
	  double cLam = pK.p3().unit().dot(pLamB.p3().unit());
	  _h_angle[1]->fill(cLam);
	  if(Lam.children().size()==2) {
	    Particle proton;
	    if(Lam.children()[0].pid()== sign*2212 &&
	       Lam.children()[1].pid()==-sign*211 ) {
	      proton = Lam.children()[0];
	    }
	    else if(Lam.children()[1].pid()== sign*2212 &&
		    Lam.children()[0].pid()==-sign*211 ){
	      proton = Lam.children()[1];
	    }
	    if(proton.pid()==sign*2212) {
	      LorentzTransform boostL =  LorentzTransform::mkFrameTransformFromBeta(pLam.betaVec());
	      FourMomentum pp = boostL.transform(boost2.transform(boost.transform(proton.momentum())));
	      double cTheta = pp.p3().unit().dot(pLam.p3().unit());
	      _h_angle[0]->fill(cTheta);
	    }
	  }
	}
	// B0 -> K*0
	else if(imode==2) {
	  const Particle & Kstar  = BB.decayProducts()[ix].at( sign*313)[0];
	  Particle KK;
	  if(Kstar.children()[0].abspid()==321 &&
	     Kstar.children()[1].abspid()==211)
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==321 &&
		  Kstar.children()[0].abspid()==211)
	    KK = Kstar.children()[1];
	  else continue;
	  FourMomentum pKstar = boost.transform(Kstar.momentum());
	  FourMomentum pK     = boost.transform(KK   .momentum());
	  const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	  pK = boost3.transform(pK);
	  FourMomentum pB = boost3.transform(boost.transform(BB.decaying()[ix].momentum()));
	  double cosK = -pB.p3().unit().dot(pK.p3().unit());
	  _h_angle[2]->fill(cosK);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix%2==0) scale(_h_mass[ix],1e6/ *_c[0]);
	else        scale(_h_mass[ix],1e6/ *_c[1]);
	normalize(_h_angle[ix]);
      }
    }

    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_angle[3];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I803343);

}
