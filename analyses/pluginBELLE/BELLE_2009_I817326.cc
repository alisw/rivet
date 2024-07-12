// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K(*) l+ l-
  class BELLE_2009_I817326 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I817326);


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
      BB.addStable(PID::K0S);
      BB.addStable(PID::K0L);
      BB.addStable( 313);
      BB.addStable(-313);
      BB.addStable( 323);
      BB.addStable(-323);
      declare(BB, "BB");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_p_FL[ix],1,1+ix,2);
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_br[ix][iy],1+ix,1+iy,1);
	  book(_p_FB[ix][iy],1+ix,1+iy,4);
	  for(unsigned int iz=0;iz<2;++iz)
	    book(_h_br_B[ix][iy][iz],
		 "TMP/h_br_"+toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),
		 refData(1+ix,1+iy,1));
	}
      };
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 130,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 321,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-321,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 310,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4CC = { { 130,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode5   = { { 323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode5CC = { {-323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode6   = { { 313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode6CC = { {-313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode7   = { { 323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode7CC = { {-323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode8   = { { 313,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode8CC = { {-313,1},{ 11,1}, {-11,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	_c[0]->fill();
	if(BB.decaying()[ix].abspid()==521) _c[1]->fill();
	else                                _c[2]->fill();
	int imode=0;
	if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) ||
	    (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)))      imode=0;
	else if (BB.modeMatches(ix,3,mode2)||BB.modeMatches(ix,3,mode2CC))    imode=1;
	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC))) imode=2;
	else if (BB.modeMatches(ix,3,mode4) || BB.modeMatches(ix,3,mode4CC))  imode=3;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode5)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode5CC))) imode=4;
	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode6)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode6CC))) imode=5;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode7)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode7CC))) imode=6;
      	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode8)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode8CC))) imode=7;
      	else continue;
	int il = imode<2 || imode==4 || imode==5 ? 13 : 11;
	const Particle & lp = BB.decayProducts()[ix].at(-il)[0];
	const Particle & lm = BB.decayProducts()[ix].at( il)[0];
	double qq = (lp.momentum()+lm.momentum()).mass2();
	for(unsigned int iy=0;iy<2;++iy) {
	  if(imode<4) {
	    _h_br[1][iy]->fill(qq);
	    if(BB.decaying()[ix].abspid()==521)
	      _h_br_B[1][iy][0]->fill(qq);
	    else
	      _h_br_B[1][iy][1]->fill(qq);
	  }
	  else {
	    _h_br[0][iy]->fill(qq);
	    if(BB.decaying()[ix].abspid()==521)
	      _h_br_B[0][iy][0]->fill(qq);
	    else
	      _h_br_B[0][iy][1]->fill(qq);
	  }
	}
	// first boost to bottom frame
	const LorentzTransform boost  =
	  LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	FourMomentum plp    = boost.transform(lp   .momentum());
	FourMomentum plm    = boost.transform(lm   .momentum());
	FourMomentum pB     = boost.transform(BB.decaying()[ix].momentum());
	// lepton stuff
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta((plp+plm).betaVec());
	plp = boost2.transform(plp);
	double cTheta = plp.p3().unit().dot(boost .transform(pB ).p3().unit());
	double AFB = cTheta>0 ? 1 : -1;
	for(unsigned int iy=0;iy<2;++iy) {
	  if(imode<4) _p_FB[1][iy]->fill(qq,AFB);
	  else        _p_FB[0][iy]->fill(qq,AFB);
	}
	// only K* for FL
	if(imode<4) continue;
	int iK = BB.decaying()[ix].abspid()==521 ? 323 : 313;
	iK *=  BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
	const Particle & Kstar = BB.decayProducts()[ix].at( iK)[0];
	FourMomentum pKstar = boost.transform(Kstar.momentum());
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
	FourMomentum pK     = boost.transform(KK   .momentum());
	const LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	pK = boost3.transform(pK);
	cTheta = pK.p3().unit().dot(boost .transform(pB ).p3().unit());
	double FL = .5*(5.*sqr(cTheta)-1.);
	for(unsigned int iy=0;iy<2;++iy) _p_FL[iy]->fill(qq,FL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // ratio of lifetimes
      double rLife = 1./1.078;
      for(unsigned int ix=0;ix<2;++ix) {
      	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h_br[ix][iy],0.5e7/ *_c[0]);
	  for(unsigned int il=0;il<2;++il) {
	    scale(_h_br_B[ix][iy][il],0.5e7/ *_c[il+1]);
	    if (il==0) scale(_h_br_B[ix][iy][il],rLife);
	  }
	  // A_I plots
	  Scatter2DPtr AI;
	  book(AI,1+ix,1+iy,3);
	  asymm(_h_br_B[ix][iy][1],_h_br_B[ix][iy][0],AI);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_br[2][2],_h_br_B[2][2][2];
    Profile1DPtr _p_FL[2],_p_FB[2][2];
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2009_I817326/d01-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2009_I817326/d01-x02-y02
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I817326);

}
