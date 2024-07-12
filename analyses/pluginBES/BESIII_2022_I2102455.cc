// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0/D+ to pions
  class BESIII_2022_I2102455 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2102455);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411 or Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");
      // histograms
      vector<unsigned int> nHist = {3,7,9,14,2,2,5,4,11};
      for(unsigned int ix=0;ix<9;++ix) {
	_h.push_back(vector<Histo1DPtr>());
	for(unsigned int iy=0;iy<nHist[ix];++iy) {
	  Histo1DPtr tmp;
	  book(tmp,ix+1,1,iy+1);
	  _h[ix].push_back(tmp);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode1   = { { 211,1}, { -211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 211,1}, { -211,1}, { 111,2} };
      static const map<PdgId,unsigned int> & mode3   = { { 211,2}, { -211,2}, { 111,1} };
      static const map<PdgId,unsigned int> & mode4   = { { 211,2}, { -211,2}, { 111,2} };
      static const map<PdgId,unsigned int> & mode5   = { { 211,2}, { -211,1}};
      static const map<PdgId,unsigned int> & mode5CC = { { 211,1}, { -211,2}};
      static const map<PdgId,unsigned int> & mode6   = { { 211,1}, {  111,2}};
      static const map<PdgId,unsigned int> & mode6CC = { {-211,1}, {  111,2}};
      static const map<PdgId,unsigned int> & mode7   = { { 211,2}, { -211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode7CC = { { 211,1}, { -211,2}, { 111,1} };
      static const map<PdgId,unsigned int> & mode8   = { { 211,3}, { -211,2}};
      static const map<PdgId,unsigned int> & mode8CC = { { 211,2}, { -211,3}};
      static const map<PdgId,unsigned int> & mode9   = { { 211,3}, { -211,2}, { 111,1} };
      static const map<PdgId,unsigned int> & mode9CC = { { 211,2}, { -211,3}, { 111,1} };
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
	// D0 -> pi+ pi- pi0
	if ( DD.modeMatches(ix,3,mode1)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // KS0 veto
	  double mpm = (pip[0].momentum()+pim[0].momentum()).mass();
	  if(mpm>.468 && mpm<.528) continue;
	  _h[0][0]->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	  _h[0][1]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	  _h[0][2]->fill(mpm);
	}
	// D0 -> pi+ pi- 2pi0 
	else if ( DD.modeMatches(ix,4,mode2)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // KS0 veto
	  FourMomentum ppm = pip[0].momentum()+pim[0].momentum(); 
	  double mpm = ppm.mass();
	  if(mpm>.468 && mpm<.528) continue;
	  FourMomentum p00 = pi0[0].momentum()+pi0[1].momentum();
	  double m00 = p00.mass();
	  if(m00>.428 && m00<.548) continue;
	  _h[1][0]->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	  _h[1][0]->fill((pip[0].momentum()+pi0[1].momentum()).mass());
	  _h[1][1]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	  _h[1][1]->fill((pim[0].momentum()+pi0[1].momentum()).mass());
	  _h[1][2]->fill(mpm);
	  _h[1][3]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	  _h[1][4]->fill((p00+pip[0].momentum()).mass());
	  _h[1][5]->fill((p00+pim[0].momentum()).mass());
	  _h[1][6]->fill((ppm+pi0[0].momentum()).mass());
	  _h[1][6]->fill((ppm+pi0[1].momentum()).mass());
	}
	// D0 -> 2pi+ 2pi- pi0 
	else if ( DD.modeMatches(ix,5,mode3)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // pi+ pi- masses and KS0 veto
	  double mpm[4];
	  bool veto=false;
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      mpm[2*ix+iy] = (pip[ix].momentum()+pim[iy].momentum()).mass(); 
	      if(mpm[2*ix+iy]>.468 && mpm[2*ix+iy]<.528) veto=true;
	    }
	  }
	  if(veto) continue;
	  // fill the histograms
	  FourMomentum ppp = pip[0].momentum()+pip[1].momentum();
	  FourMomentum pmm = pim[0].momentum()+pim[1].momentum();
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h[2][0]->fill((pip[ix].momentum()+pi0[0].momentum()).mass());
	    _h[2][1]->fill((pim[ix].momentum()+pi0[0].momentum()).mass());
	    _h[2][3]->fill((ppp+pim[ix].momentum()).mass());
	    _h[2][4]->fill((pmm+pip[ix].momentum()).mass());
	    _h[2][6]->fill((ppp+pim[ix].momentum()+pi0[0].momentum()).mass());
	    _h[2][7]->fill((pmm+pip[ix].momentum()+pi0[0].momentum()).mass());
	    for(unsigned int iy=0;iy<2;++iy) {
	      _h[2][2]->fill(mpm[2*ix+iy]);
	      _h[2][5]->fill((pip[ix].momentum()+pim[iy].momentum()+pi0[0].momentum()).mass()); 
	    }
	  }
	  _h[2][8]->fill((ppp+pmm).mass());
	}
	// D0 -> 2pi+ 2pi- pi0 
	else if ( DD.modeMatches(ix,6,mode4)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // pi+ pi- masses and KS0 veto
	  double mpm[4];
	  bool veto=false;
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      mpm[2*ix+iy] = (pip[ix].momentum()+pim[iy].momentum()).mass(); 
	      if(mpm[2*ix+iy]>.468 && mpm[2*ix+iy]<.528) veto=true;
	    }
	  }
	  if(veto) continue;
	  FourMomentum p00 = pi0[0].momentum()+pi0[1].momentum();
	  double m00 = p00.mass();
	  if(m00>.428 && m00<.548) continue;
	  // fill the histograms
	  _h[3][3]->fill(m00);
	  FourMomentum ppp = pip[0].momentum()+pip[1].momentum();
	  FourMomentum pmm = pim[0].momentum()+pim[1].momentum();
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int i0=0;i0<2;++i0) {
	      _h[3][ 0]->fill((pip[ix].momentum()+pi0[i0].momentum()).mass());
	      _h[3][ 1]->fill((pim[ix].momentum()+pi0[i0].momentum()).mass());
	      _h[3][ 7]->fill((ppp+pim[ix].momentum()+pi0[i0].momentum()).mass());
	      _h[3][ 8]->fill((pmm+pip[ix].momentum()+pi0[i0].momentum()).mass());
	      _h[3][11]->fill((ppp+pim[ix].momentum()+p00).mass());
	      _h[3][12]->fill((pmm+pip[ix].momentum()+p00).mass());
	      _h[3][13]->fill((pmm+pi0[ix].momentum()+ppp).mass());
	      for(unsigned int iy=0;iy<2;++iy)
		_h[3][6]->fill((pip[ix].momentum()+pim[iy].momentum()+pi0[i0].momentum()).mass());
	    }
	    _h[3][4]->fill((ppp+pim[ix].momentum()).mass());
	    _h[3][5]->fill((pmm+pip[ix].momentum()).mass());
	    for(unsigned int iy=0;iy<2;++iy) {
	      _h[3][ 2]->fill(mpm[2*ix+iy]);
	      _h[3][10]->fill((pip[ix].momentum()+pim[iy].momentum()+p00).mass());
	    }
	  }
	  _h[3][9]->fill((ppp+pmm).mass());
	}
	// D+ -> 2pi+ pi- 
	else if ( DD.modeMatches(ix,3,mode5  ) ||
		  DD.modeMatches(ix,3,mode5CC)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  // pi+ pi- masses and KS0 veto
	  double mpm[2];
	  bool veto=false;
	  for(unsigned int ix=0;ix<2;++ix) {
	    mpm[ix] = (pip[ix].momentum()+pim[0].momentum()).mass(); 
	    if(mpm[ix]>.468 && mpm[ix]<.528) veto=true;
	  }
	  if(veto) continue;
	  _h[4][0]->fill((pip[0].momentum()+pip[1].momentum()).mass());
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h[4][1]->fill(mpm[ix]);
	  }
	}
	// D+ -> pi+ 2pi0 
	else if ( DD.modeMatches(ix,3,mode6  ) ||
		  DD.modeMatches(ix,3,mode6CC)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  FourMomentum p00 = pi0[0].momentum()+pi0[1].momentum();
	  double m00 = p00.mass();
	  if(m00>.428 && m00<.548) continue;
	  // fill the histograms
	  _h[5][1]->fill(m00);
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h[5][0]->fill((pip[0].momentum()+pi0[ix].momentum()).mass());
	  }
	}
	// D+ -> 2pi+ pi- pi0
	else if ( DD.modeMatches(ix,4,mode7  ) ||
		  DD.modeMatches(ix,4,mode7CC)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // pi+ pi- masses and KS0 veto
	  double mpm[2];
	  bool veto=false;
	  for(unsigned int ix=0;ix<2;++ix) {
	    mpm[ix] = (pip[ix].momentum()+pim[0].momentum()).mass(); 
	    if(mpm[ix]>.468 && mpm[ix]<.528) veto=true;
	  }
	  if(veto) continue;
	  _h[6][1]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	  _h[6][3]->fill((pim[0].momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h[6][0]->fill((pip[ix].momentum()+pi0[0].momentum()).mass());
	    _h[6][2]->fill(mpm[ix]);
	    _h[6][4]->fill((pip[ix].momentum()+pim[0].momentum()+pi0[0].momentum()).mass());
	  }
	}
	// D+ -> 3pi+ 2pi- 
	else if ( DD.modeMatches(ix,5,mode8  ) ||
		  DD.modeMatches(ix,5,mode8CC)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  // pi+ pi- masses and KS0 veto
	  double mpm[3][2];
	  bool veto=false;
	  for(unsigned int ix=0;ix<3;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      mpm[ix][iy] = (pip[ix].momentum()+pim[iy].momentum()).mass(); 
	      if(mpm[ix][iy]>.468 && mpm[ix][iy]<.528) veto=true;
	    }
	  }
	  if(veto) continue;
	  FourMomentum pppp = pip[0].momentum()+pip[1].momentum()+pip[2].momentum();
	  FourMomentum pmm  = pim[0].momentum()+pim[1].momentum();
	  for(unsigned int ix=0;ix<3;++ix) {
	    _h[7][2]->fill((pmm+pip[ix].momentum()).mass());
	    for(unsigned int iy=0;iy<2;++iy) {
	      _h[7][0]->fill(mpm[ix][iy]);
	      _h[7][1]->fill((pppp-pip[ix].momentum()+pim[iy].momentum()).mass());
	    }
	  }
	  for(unsigned int iy=0;iy<2;++iy)
	    _h[7][3]->fill((pppp+pim[iy].momentum()).mass());
	}
	// D+ -> 3pi+ 2pi- pi0
	else if ( DD.modeMatches(ix,6,mode9  ) ||
		  DD.modeMatches(ix,6,mode9CC)) {
	  const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	  const Particles & pim = DD.decayProducts()[ix].at(-sign*211);
	  const Particles & pi0 = DD.decayProducts()[ix].at(      111);
	  // pi+ pi- masses and KS0 veto
	  double mpm[3][2];
	  bool veto=false;
	  for(unsigned int ix=0;ix<3;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      mpm[ix][iy] = (pip[ix].momentum()+pim[iy].momentum()).mass(); 
	      if(mpm[ix][iy]>.468 && mpm[ix][iy]<.528) veto=true;
	    }
	  }
	  if(veto) continue;
	  FourMomentum pppp = pip[0].momentum()+pip[1].momentum()+pip[2].momentum();
	  FourMomentum pmm  = pim[0].momentum()+pim[1].momentum();
	  FourMomentum pcharged = pppp+pmm;
	  _h[8][9]->fill(pcharged.mass());
	  for(unsigned int ix=0;ix<3;++ix) {
	    _h[8][0]->fill((pip[ix].momentum()+pi0[0].momentum()).mass());
	    FourMomentum pmmp = (pmm+pip[ix].momentum()); 
	    _h[8][ 4]->fill(pmmp.mass()); 
	    _h[8][ 7]->fill((pmmp+pi0[0].momentum()).mass());
	    _h[8][ 8]->fill((pcharged-pip[ix].momentum()).mass());
	    _h[8][10]->fill((pcharged-pip[ix].momentum()+pi0[0].momentum()).mass());
	    for(unsigned int iy=0;iy<2;++iy) {
	      _h[8][2]->fill(mpm[ix][iy]);
	      FourMomentum pppm = (pppp-pip[ix].momentum()+pim[iy].momentum()); 
	      _h[8][3]->fill(pppm.mass()); 
	      _h[8][6]->fill((pppm+pi0[0].momentum()).mass());
	      _h[8][5]->fill((pip[ix].momentum()+pim[iy].momentum()+pi0[0].momentum()).mass());
	    }
	  }
	  for(unsigned int iy=0;iy<2;++iy) {
	    _h[8][1]->fill((pim[iy].momentum()+pi0[0].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<_h.size();++ix) {
	for(unsigned int iy=0;iy<_h[ix].size();++iy) {
	  normalize(_h[ix][iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    vector<vector<Histo1DPtr> > _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2102455);

}
