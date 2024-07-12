// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> p pbar + pions
  class BABAR_2012_I946659 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I946659);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 411);
      BB.addStable(-411);
      BB.addStable( 421);
      BB.addStable(-421);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( 423);
      BB.addStable(-423);
      declare(BB, "BB");
      // histograms
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int ix=0;ix<4;++ix) {
	  if(ix<2) book(_h[ix][iy],1,1+ix,1+iy);
	  book(_h[ix+2][iy],2,1+ix,1+iy);
	  book(_h[ix+6][iy],3,1+ix,1+iy);
	}
      }
      book(_nB[0],"/TMP/nB0");
      book(_nB[1],"/TMP/nBP");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 421,1},{ 2212,1}, {-2212,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-421,1},{ 2212,1}, {-2212,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 423,1},{ 2212,1}, {-2212,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-423,1},{ 2212,1}, {-2212,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 411,1},{ 2212,1}, {-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-411,1},{ 2212,1}, {-2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 413,1},{ 2212,1}, {-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode4CC = { {-413,1},{ 2212,1}, {-2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode5   = { { 421,1},{ 2212,1}, {-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode5CC = { {-421,1},{ 2212,1}, {-2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode6   = { { 423,1},{ 2212,1}, {-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode6CC = { {-423,1},{ 2212,1}, {-2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode7    = { { 421,1},{ 2212,1}, {-2212,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode7CC  = { {-421,1},{ 2212,1}, {-2212,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode8    = { { 423,1},{ 2212,1}, {-2212,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode8CC  = { {-423,1},{ 2212,1}, {-2212,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode9    = { { 411,1},{ 2212,1}, {-2212,1}, { 211,2}};
      static const map<PdgId,unsigned int> & mode9CC  = { {-411,1},{ 2212,1}, {-2212,1}, { 211,2}};
      static const map<PdgId,unsigned int> & mode10   = { { 413,1},{ 2212,1}, {-2212,1}, {-211,2}};
      static const map<PdgId,unsigned int> & mode10CC = { {-413,1},{ 2212,1}, {-2212,1}, { 211,2}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	int sign = 1, imode = -1, idd=421;
	if(BB.decaying()[ix].abspid()==511) {
	  _nB[0]->fill();
	  if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1)) {
	    sign=1;
	    imode=0;
	    idd=421;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1CC)) {
	    sign=-1;
	    imode=0;
	    idd=421;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2)) {
	    sign=1;
	    imode=1;
	    idd=423;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2CC)) {
	    sign=-1;
	    imode=1;
	    idd=423;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode3)) {
	    sign=1;
	    imode=2;
	    idd=411;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode3CC)) {
	    sign=-1;
	    imode=2;
	    idd=411;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode4)) {
	    sign=1;
	    imode=3;
	    idd=413;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode4CC)) {
	    sign=-1;
	    imode=3;
	    idd=413;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,5,mode7)) {
	    sign=1;
	    imode=6;
	    idd=421;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,5,mode7CC)) {
	    sign=-1;
	    imode=6;
	    idd=421;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,5,mode8)) {
	    sign=1;
	    imode=7;
	    idd=423;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,5,mode8CC)) {
	    sign=-1;
	    imode=7;
	    idd=423;
	  }
	  else
	    continue;
	}
	else {
	  _nB[1]->fill();
	  if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode5)) {
	    sign=1;
	    imode=4;
	    idd=421;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode5CC)) {
	    sign=-1;
	    imode=4;
	    idd=421;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode6)) {
	    sign=1;
	    imode=5;
	    idd=423;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode6CC)) {
	    sign=-1;
	    imode=5;
	    idd=423;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,5,mode9)) {
	    sign=1;
	    imode=8;
	    idd=411;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,5,mode9CC)) {
	    sign=-1;
	    imode=8;
	    idd=411;
	  }
	  else if (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,5,mode10)) {
	    sign=1;
	    imode=9;
	    idd=413;
	  }
	  else if  (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,5,mode10CC)) {
	    sign=-1;
	    imode=9;
	    idd=413;
	  }
	  else
	    continue;
	}
	const Particle & DD   = BB.decayProducts()[ix].at( sign*idd )[0];
	const Particle & pp   = BB.decayProducts()[ix].at( sign*2212)[0];
	const Particle & pbar = BB.decayProducts()[ix].at(-sign*2212)[0];
	double mppbar = (pp.momentum()+pbar.momentum()).mass();
	double mpD    = (pp.momentum()+DD.momentum()).mass();
	if(imode==0) {
	  if(sqr(mppbar)>5.) _h[0][0]->fill(mpD);
	  else               _h[0][1]->fill(mpD);
	  if(sqr(mpD)>9.)    _h[0][2]->fill(mppbar);
	  else               _h[0][3]->fill(mppbar);
	}
	else if(imode==1) {
	  if(sqr(mppbar)>5.) _h[1][0]->fill(mpD);
	  else               _h[1][1]->fill(mpD);
	  if(sqr(mpD)>10.5)  _h[1][2]->fill(mppbar);
	  else               _h[1][3]->fill(mppbar);
	}
	else {
	  _h[imode][0]->fill(mppbar);
	  _h[imode][1]->fill((pbar.momentum()+DD.momentum()).mass());
	  _h[imode][2]->fill(mpD);
	  const Particles & pim   = BB.decayProducts()[ix].at(-sign*211);
	  for(unsigned int iy=0;iy<pim.size();++iy) {
	    _h[imode][3]->fill((pp.momentum()+pim[iy].momentum()).mass());
	  }
	  if(imode==6 || imode==7 ) {
	    const Particle & pip = BB.decayProducts()[ix].at(sign*211)[0];
	    _h[imode][3]->fill((pp.momentum()+pip.momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<10;++ix) {
	CounterPtr temp = (ix<4 || ix==6 || ix==7) ? _nB[0] : _nB[1];
	for(unsigned int iy=0;iy<4;++iy) {
	  if(iy==3 && ix>5)
	    scale(_h[ix][iy],0.5e5/ *temp);
	  else
	    scale(_h[ix][iy],1e5/ *temp);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[10][4];
    CounterPtr _nB[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I946659);

}
