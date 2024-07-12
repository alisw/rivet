// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K K K decays
  class BABAR_2012_I1086537 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1086537);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h_aver[ix][iy],1+2*ix,1,1+iy);
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_charge2[ix],4,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_charge1[ix][iy],2,1+iy,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,2},{-321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 321,1},{-321,2}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,2},{ 321,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 310,2},{-321,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 310,1},{ 321,1} ,{-321,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int sign = 1, imode = 0;
      	if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	  imode=0;
      	  sign=1;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
	  imode=0;
      	  sign=-1;
      	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	  imode=1;
      	  sign=1;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
	  imode=1;
      	  sign=-1;
      	}
	else if (BB.modeMatches(ix,3,mode3)) {
	  imode=2;
	  if(BB.decaying()[ix].pid()<0) sign=-1;
	}
      	else
      	  continue;
	// B+ -> K+ K+ K-
	if(imode==0) {
	  const Particles & Kp = BB.decayProducts()[ix].at( sign*321);
	  const Particle  & Km = BB.decayProducts()[ix].at(-sign*321)[0];
	  double mKpKm[2];
	  for(unsigned int ix=0;ix<2;++ix)
	    mKpKm[ix] = (Kp[ix].momentum()+Km.momentum()).mass();
	  if(mKpKm[0]>mKpKm[1]) swap(mKpKm[0],mKpKm[1]);
	  _h_aver[0][0]->fill(mKpKm[0]);
	  _h_aver[0][1]->fill(mKpKm[1]);
	  _h_aver[0][2]->fill((Kp[0].momentum()+Kp[1].momentum()).mass());
	  if(BB.decaying()[ix].pid()>0) {
	    for(unsigned int ix=0;ix<2;++ix) _h_charge1[0][ix]->fill(mKpKm[0]);
	  }
	  else {
	    for(unsigned int ix=0;ix<2;++ix) _h_charge1[1][ix]->fill(mKpKm[0]);
	  }
	}
	// B+ -> KS0 KS0 K+
	else if(imode==1) {
	  const Particles & K0 = BB.decayProducts()[ix].at( 310);
	  const Particle  & Kp = BB.decayProducts()[ix].at( sign*321)[0];
	  double mK0K0 = (K0[0].momentum()+K0[1].momentum()).mass();
	  _h_aver[1][0]->fill(mK0K0);
	  double mKpK0[2];
	  for(unsigned int ix=0;ix<2;++ix)
	    mKpK0[ix] = (K0[ix].momentum()+Kp.momentum()).mass();
	  if(mKpK0[0]>mKpK0[1]) swap(mKpK0[0],mKpK0[1]);
	  _h_aver[1][1]->fill(mKpK0[1]);
	  _h_aver[1][2]->fill(mKpK0[0]);
	  if(BB.decaying()[ix].pid()>0) {
	    _h_charge2[0]->fill(mK0K0);
	  }
	  else {
	    _h_charge2[1]->fill(mK0K0);
	  }
	}
	// B0 -> KS0 K+ K-
	else if(imode==2) {
	  const Particle & K0 = BB.decayProducts()[ix].at(      310)[0];
	  const Particle & Kp = BB.decayProducts()[ix].at( sign*321)[0];
	  const Particle & Km = BB.decayProducts()[ix].at(-sign*321)[0];
	  _h_aver[2][0]->fill((Kp.momentum()+Km.momentum()).mass());
	  _h_aver[2][1]->fill((K0.momentum()+Km.momentum()).mass());
	  _h_aver[2][2]->fill((K0.momentum()+Kp.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h_aver[ix][iy],1.,false);
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_charge2[ix],1.,false);
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h_charge1[ix][iy],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_aver[3][3];
    Histo1DPtr _h_charge1[2][2];
    Histo1DPtr _h_charge2[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1086537);

}
