// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_c - > jpsi + 3 charged hadrons
  class LHCB_2022_I1960979 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2022_I1960979);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==541);
      declare(ufs, "UFS");
      DecayedParticles BC(ufs);
      BC.addStable( PID::PI0);
      BC.addStable( PID::K0S);
      BC.addStable( PID::JPSI);
      declare(BC, "BC");
      for(unsigned int ix=0;ix<4;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  if(iy==1&&ix>1) continue;
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 443,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 443,1}, {-211,2}, { 211,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 443,1}, { 321,1}, {-321,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode2CC = { { 443,1}, { 321,1}, {-321,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode3   = { { 443,1}, { 211,1}, {-211,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode3CC = { { 443,1}, { 211,1}, {-211,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode4   = { { 443,1}, { 321,2}, {-321,1} };
      static const map<PdgId,unsigned int> & mode4CC = { { 443,1}, { 321,1}, {-321,2} };
      DecayedParticles BC = apply<DecayedParticles>(event, "BC");
      // loop over particles
      for(unsigned int ix=0;ix<BC.decaying().size();++ix) {
	int sign = BC.decaying()[ix].pid()/BC.decaying()[ix].abspid();
	if ((sign== 1 && BC.modeMatches(ix,4,mode1  )) ||
	    (sign==-1 && BC.modeMatches(ix,4,mode1CC))) {
	  const Particle  & pim = BC.decayProducts()[ix].at(-sign*211)[0];
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  _h[0][0]->fill((pim.momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	  _h[0][1]->fill((pim.momentum()+pip[0].momentum()).mass());
	  _h[0][1]->fill((pim.momentum()+pip[1].momentum()).mass());
	}
	else if((sign== 1 && BC.modeMatches(ix,4,mode2  )) ||
		(sign==-1 && BC.modeMatches(ix,4,mode2CC))) {
	  const Particle  & Kp  = BC.decayProducts()[ix].at( sign*321)[0];
	  const Particle  & Km  = BC.decayProducts()[ix].at(-sign*321)[0];
	  const Particle  & pip = BC.decayProducts()[ix].at( sign*211)[0];
	  _h[1][0]->fill((Km.momentum()+pip.momentum()).mass());
	  _h[1][1]->fill((Km.momentum()+Kp .momentum()).mass());
	}
	else if((sign== 1 && BC.modeMatches(ix,4,mode3  )) ||
		(sign==-1 && BC.modeMatches(ix,4,mode3CC))) {
	  const Particle  & pim  = BC.decayProducts()[ix].at(-sign*211)[0];
	  const Particle  & Kp   = BC.decayProducts()[ix].at( sign*321)[0];
	  _h[2][0]->fill((Kp.momentum()+pim.momentum()).mass());
	}
	else if((sign== 1 && BC.modeMatches(ix,4,mode4  )) ||
		(sign==-1 && BC.modeMatches(ix,4,mode4CC))) {
	  const Particle  & Km = BC.decayProducts()[ix].at(-sign*321)[0];
	  const Particles & Kp = BC.decayProducts()[ix].at( sign*321);
	  _h[3][0]->fill((Kp[0].momentum()+Km.momentum()).mass());
	  _h[3][0]->fill((Kp[1].momentum()+Km.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  if(iy==1&&ix>1) continue;
	  normalize(_h[ix][iy],1.,false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2022_I1960979);

}
