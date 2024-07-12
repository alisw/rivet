// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K(*) l+ l-
  class BABAR_2012_I1111233 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1111233);


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
      BB.addStable( 310);
      BB.addStable( 313);
      BB.addStable( 323);
      BB.addStable(-313);
      BB.addStable(-323);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_br[ix][iy],1,1+ix,1+iy);
	  book(_p_CP[ix][iy],2,1+ix,1+iy);
	  for(unsigned int il=0;il<2;++il) {
	    book(_h_br_l[il][ix][iy],"TMP/h_br_l_"+toString(il)+"_"+toString(ix)+"_"+toString(iy),
		 refData(3,1+ix,1+iy));
	    book(_h_br_I[il][ix][iy],"TMP/h_br_I_"+toString(il)+"_"+toString(ix)+"_"+toString(iy),
		 refData(4,1+ix,1+iy));
	  }
	}
      }
      for(unsigned int ix=0;ix<3;++ix)
	book(_c[ix],"TMP/c_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // kaon modes
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 321,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-321,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 310,1},{ 11,1}, {-11,1}};
      // K* modes
      static const map<PdgId,unsigned int> & mode5   = { { 323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode5CC = { {-323,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode6   = { { 313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode6CC = { {-313,1},{ 13,1}, {-13,1}};
      static const map<PdgId,unsigned int> & mode7   = { { 323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode7CC = { {-323,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode8   = { { 313,1},{ 11,1}, {-11,1}};
      static const map<PdgId,unsigned int> & mode8CC = { {-313,1},{ 11,1}, {-11,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	_c[0]->fill();
	if(BB.decaying()[ix].abspid()==521) _c[1]->fill();
	else                                _c[2]->fill();
	int imode=0;
	if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) ||
	    (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)))      imode=0;
	else if (BB.modeMatches(ix,3,mode2))                                  imode=1;
	else if ((BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) ||
		 (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC))) imode=2;
	else if (BB.modeMatches(ix,3,mode4))                                  imode=3;
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
	double ACP = BB.decaying()[ix].pid()>0 ? -1. : 1.;
	if(imode<4) {
	  double wgt = (imode==1||imode==3) ? 2 : 1;
	  for(unsigned int iy=0;iy<2;++iy) {
	    _h_br[iy][0]->fill(qq,wgt);
	    _p_CP[iy][0]->fill(qq,ACP,wgt);
	    if(il==13) _h_br_l[0][iy][0]->fill(qq,wgt);
	    else       _h_br_l[1][iy][0]->fill(qq,wgt);
	    if(BB.decaying()[ix].abspid()==521)
	      _h_br_I[0][iy][0]->fill(qq,wgt);
	    else
	      _h_br_I[1][iy][0]->fill(qq,wgt);
	  }
	}
	else {
	  for(unsigned int iy=0;iy<2;++iy) {
	    _h_br[iy][1]->fill(qq);
	    _p_CP[iy][1]->fill(qq,ACP);
	    if(BB.decaying()[ix].abspid()==521)
	      _h_br_I[0][iy][1]->fill(qq);
	    else
	      _h_br_I[1][iy][1]->fill(qq);
	    if(il==13) _h_br_l[0][iy][1]->fill(qq);
	    else       _h_br_l[1][iy][1]->fill(qq);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // ratio of lifetimes
      double rLife = 1./1.078;
      for(unsigned int ix=0;ix<2;++ix) {
      	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h_br[ix][iy],1e7/ *_c[0]);
	  for(unsigned int il=0;il<2;++il) {
	    scale(_h_br_l[il][ix][iy],1e7/ *_c[0]);
	    scale(_h_br_I[il][ix][iy],1e7/ *_c[il+1]);
	    if (il==0) scale(_h_br_I[il][ix][iy],rLife);
	  }
	  // RK plots
	  Scatter2DPtr RK;
	  book(RK,3,1+ix,1+iy);
	  divide(_h_br_l[0][ix][iy],_h_br_l[1][ix][iy],RK);
	  book(RK,4,1+ix,1+iy);
	  asymm(_h_br_I[1][ix][iy],_h_br_I[0][ix][iy],RK);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_br[2][2],_h_br_l[2][2][2],_h_br_I[2][2][2];
    Profile1DPtr _p_CP[2][2];
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1111233);

}
