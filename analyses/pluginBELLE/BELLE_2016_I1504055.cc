// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B -> K* l+l-
  class BELLE_2016_I1504055 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1504055);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(   443);
      BB.addStable(100443);
      BB.addStable( 313);
      BB.addStable( 323);
      BB.addStable(-313);
      BB.addStable(-323);
      declare(BB, "BB");
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<6;++iy) {
	  book(_p_P[ix][iy],1,1+ix,1+iy);
	  if(iy>1) continue;
	  book(_p_Q[ix][iy],2,1+ix,1+iy);
	}
      }
      book(_FL,"TMP/FL");
      book(_norm,"TMP/norm");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 313,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4CC = { {-313,1},{ 11,1}, {-11,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int imode=0;
      	if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) ||
	    (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)))       imode=0;
	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)))  imode=1;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC)))  imode=2;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode4)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode4CC)))  imode=3;
      	else continue;
      	int il = imode<2 ? 13 : 11;
	int sign = BB.decaying()[ix].pid()>0 ? 1 : -1;
      	const Particle & lp = BB.decayProducts()[ix].at(-sign*il)[0];
      	const Particle & lm = BB.decayProducts()[ix].at( sign*il)[0];
      	double qq = (lp.momentum()+lm.momentum()).mass2();
	int iK = BB.decaying()[ix].abspid()==521 ? 323 : 313;
	iK *=  BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
	const Particle & Kstar = BB.decayProducts()[ix].at( iK)[0];
	if(Kstar.children().size()!=2) continue;
	Particle KK;
	if(Kstar.abspid()==313) {
	  if(Kstar.children()[0].abspid()==321 &&
	     Kstar.children()[1].abspid()==211)
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==321 &&
		  Kstar.children()[0].abspid()==211)
	    KK = Kstar.children()[1];
	  else continue;
	}
	else {
	  if(Kstar.children()[0].abspid()==311 &&
	     Kstar.children()[1].abspid()==211)
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==311 &&
		  Kstar.children()[0].abspid()==211)
	    KK = Kstar.children()[1];
	  else if(Kstar.children()[0].abspid()==310 &&
		  Kstar.children()[1].abspid()==211)
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==310 &&
		  Kstar.children()[0].abspid()==211)
	    KK = Kstar.children()[1];
	  else if(Kstar.children()[0].abspid()==321 &&
		  Kstar.children()[1].abspid()==111 )
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==321 &&
		  Kstar.children()[0].abspid()==111 )
	    KK = Kstar.children()[1];
	  else continue;
	  if(KK.abspid()==311) {
	    if(KK.children().size()==1 && KK.children()[0].pid()==310)
	      KK = KK.children()[0];
	    else
	      continue;
	  }
	}
	// first boost to bottom frame
	const LorentzTransform boost  = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	FourMomentum plp    = boost.transform(lp   .momentum());
	FourMomentum plm    = boost.transform(lm   .momentum());
	FourMomentum pKstar = boost.transform(Kstar.momentum());
	FourMomentum pK     = boost.transform(KK   .momentum());
	FourMomentum pB     = boost.transform(BB.decaying()[ix].momentum());
	// lepton stuff
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta((plp+plm).betaVec());
	plp = boost2.transform(plp);
	Vector3 axis1 = boost .transform(pB ).p3().unit();
	double cThetaL = plp.p3().unit().dot(axis1);
	Vector3 Trans1 = plp.p3() - cThetaL*plp.p3().mod()*axis1;
	// kaon stuff
	const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	pK = boost3.transform(pK);
	Vector3 axis2 = boost .transform(pB ).p3().unit();
	double cThetaK = pK.p3().unit().dot(axis2);
	double FL = .5*(5.*sqr(cThetaK)-1.);
	Vector3 Trans2 = pK.p3() - cThetaK*pK.p3().mod()*axis2;
	double phi = atan2(Trans1.cross(Trans2).dot(axis2),Trans1.dot(Trans2));
	double sThetaL = sqrt(1.-sqr(cThetaL));
	double sThetaK = sqrt(1.-sqr(cThetaK));
	double S4 = 12.5*cThetaL*sThetaL*cThetaK*sThetaK*cos(phi);
	double S5 = 5.*cThetaK*sThetaK*sThetaL*sin(phi);
	_FL->fill(FL);
	_norm->fill();
	for(unsigned int ix=0;ix<2;++ix) {
	  _p_P[ix][0]->fill(qq,S4);
	  _p_P[ix][3]->fill(qq,S5);
	  if(il==11) {
	    _p_P[ix][1]->fill(qq,S4);
	    _p_P[ix][4]->fill(qq,S5);
	    _p_Q[ix][0]->fill(qq,-S4);
	    _p_Q[ix][1]->fill(qq,-S5);
	  }
	  else {
	    _p_P[ix][2]->fill(qq,S4);
	    _p_P[ix][5]->fill(qq,S5);
	    _p_Q[ix][0]->fill(qq,S4);
	    _p_Q[ix][1]->fill(qq,S5);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter1D R = *_FL/ *_norm;
      double fl = R.point(0).x();
      double fact = 1./sqrt(fl*(1.-fl));
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<6;++iy) {
	  _p_P[ix][iy]->scaleY(fact);
	  if(iy>1) continue;
	  _p_Q[ix][iy]->scaleY(fact);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p_P[2][6],_p_Q[2][2];
    CounterPtr _FL,_norm;
    /// @}
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y03
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y04
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y05
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x01-y06
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y03
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y04
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y05
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d01-x02-y06
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d02-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d02-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d02-x02-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1504055/d02-x02-y02


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1504055);

}
