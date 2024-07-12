// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K* l+l-
  class BABAR_2016_I1391152 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2016_I1391152);


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
      // book histograms
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_p_FL[ix][iy],1,1+ix,1+iy);
	  book(_p_FB[ix][iy],2,1+ix,1+iy);
	  book(_p_P2_num[ix][iy],"TMP/P2_num_"+toString(ix+1)+"_"+toString(iy+1),refData(3,1+ix,1+iy));
	  book(_p_P2_den[ix][iy],"TMP/P2_den_"+toString(ix+1)+"_"+toString(iy+1),refData(3,1+ix,1+iy));
	}
      }
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
		  Kstar.children()[1].abspid()==111 && il==11)
	    KK = Kstar.children()[0];
	  else if(Kstar.children()[1].abspid()==321 &&
		  Kstar.children()[0].abspid()==111 && il==11)
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
	double cTheta = plp.p3().unit().dot(boost .transform(pB ).p3().unit());
	double AFB = cTheta>0 ? 1 : -1;
	// kaon stuff
	const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	pK = boost3.transform(pK);
	cTheta = pK.p3().unit().dot(boost .transform(pB ).p3().unit());
	double FL = .5*(5.*sqr(cTheta)-1.);
	// fill histograms
	for(unsigned int iz=0;iz<2;++iz) {
	  for(unsigned int iy=0;iy<3;++iy) {
	    if( (BB.decaying()[ix].abspid()==511&&iy==0) ||
		(BB.decaying()[ix].abspid()==521&&iy==1) ) continue;
	    _p_FB    [iz][iy]->fill(qq, AFB);
	    _p_FL    [iz][iy]->fill(qq,FL);
	    _p_P2_num[iz][iy]->fill(qq, -2./3.*AFB);
	    _p_P2_den[iz][iy]->fill(qq, 1.-FL);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  Scatter2DPtr tmp;
	  book(tmp,3,1+ix,1+iy);
	  divide(_p_P2_num[ix][iy],_p_P2_den[ix][iy],tmp);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p_FL[2][3],_p_FB[2][3],_p_P2_num[2][3],_p_P2_den[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2016_I1391152);

}
