// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_c to ccbar + hadrons
  class LHCB_2022_I2138845 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2022_I2138845);


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
      BC.addStable( PID::PSI2S);
      declare(BC, "BC");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix+3],2,1,1+ix);
	book(_h[ix+5],3,1,1+ix);
	book(_h[ix+7],4,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 443,1}, { 211,3}, {-211,2} };
      static const map<PdgId,unsigned int> & mode1CC = { { 443,1}, {-211,3}, { 211,2} };
      static const map<PdgId,unsigned int> & mode2   = { { 443,1}, { 321,1}, {-321,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2CC = { { 443,1}, { 321,1}, {-321,1}, {-211,2}, { 211,1} };
      static const map<PdgId,unsigned int> & mode3   = { { 443,1}, { 211,4}, {-211,3} };
      static const map<PdgId,unsigned int> & mode3CC = { { 443,1}, {-211,4}, { 211,3} };
      static const map<PdgId,unsigned int> & mode4   = { { 100443,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode4CC = { { 100443,1}, {-211,2}, { 211,1} };
      DecayedParticles BC = apply<DecayedParticles>(event, "BC");
      // loop over particles
      for(unsigned int ix=0;ix<BC.decaying().size();++ix) {
	int sign = BC.decaying()[ix].pid()/BC.decaying()[ix].abspid();
	if ((sign== 1 && BC.modeMatches(ix,6,mode1  )) ||
	    (sign==-1 && BC.modeMatches(ix,6,mode1CC))) {
	  const Particles & pim = BC.decayProducts()[ix].at(-sign*211);
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  FourMomentum ptotal=pim[0].momentum()+pim[1].momentum()+
	    pip[0].momentum()+pip[1].momentum()+pip[2].momentum();
	  _h[0]->fill(ptotal.mass());
	  for(unsigned int ix=0;ix<3;++ix)
	    for(unsigned int iy=0;iy<2;++iy) {
	      FourMomentum prho=pip[ix].momentum()+pim[iy].momentum();
	      _h[3]->fill((ptotal-prho).mass());
	      _h[5]->fill(prho.mass());
	    }
	}
	else if ((sign== 1 && BC.modeMatches(ix,6,mode2  )) ||
		 (sign==-1 && BC.modeMatches(ix,6,mode2CC))) {
	  const Particles & pim = BC.decayProducts()[ix].at(-sign*211);
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  const Particle  & Kp  = BC.decayProducts()[ix].at( sign*321)[0];
	  const Particle  & Km  = BC.decayProducts()[ix].at(-sign*321)[0];
	  _h[1]->fill((Kp.momentum()+Km.momentum()+pim[0].momentum()+
		       pip[0].momentum()+pip[1].momentum()).mass());
	  _h[7]->fill((Kp.momentum()+pim[0].momentum()).mass());
	  _h[7]->fill((Km.momentum()+pip[0].momentum()).mass());
	  _h[7]->fill((Km.momentum()+pip[1].momentum()).mass());
	  _h[8]->fill((Kp.momentum()+Km.momentum()).mass());
	}
	else if ((sign== 1 && BC.modeMatches(ix,8,mode3  )) ||
		 (sign==-1 && BC.modeMatches(ix,8,mode3CC))) {
	  const Particles & pim = BC.decayProducts()[ix].at(-sign*211);
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  _h[2]->fill((pim[0].momentum()+pim[1].momentum()+pim[2].momentum()+
		       pip[0].momentum()+pip[1].momentum()+pip[2].momentum()+pip[3].momentum()).mass());
	}
	else if ((sign== 1 && BC.modeMatches(ix,4,mode4  )) ||
		 (sign==-1 && BC.modeMatches(ix,4,mode4CC))) {
	  const Particles & pim = BC.decayProducts()[ix].at(-sign*211);
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  _h[4]->fill((pim[0].momentum()+
		       pip[0].momentum()+pip[1].momentum()).mass());
	  for(unsigned int iy=0;iy<2;++iy)
	    _h[6]->fill((pim[0].momentum()+pip[ix].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2022_I2138845);

}
