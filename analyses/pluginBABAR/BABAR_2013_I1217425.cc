// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> Lambda_c+ pbar pi+ pi-
  class BABAR_2013_I1217425 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2013_I1217425);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      // treat Lambda_c, Sigma_c and sigma_c* as stable
      B0.addStable( 4122);
      B0.addStable(-4122);
      B0.addStable( 4112);
      B0.addStable(-4112);
      B0.addStable( 4212);
      B0.addStable(-4212);
      B0.addStable( 4222);
      B0.addStable(-4222);
      B0.addStable( 4114);
      B0.addStable(-4114);
      B0.addStable( 4214);
      B0.addStable(-4214);
      B0.addStable( 4224);
      B0.addStable(-4224);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<5;++ix) {
	for(unsigned int iy=0;iy<6;++iy) {
	  if( (ix<3 && iy>=3) or (ix==4 && iy>=4) ) continue;
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }

    
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 4112,1},{-2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-4112,1},{ 2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 4222,1},{-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-4222,1},{ 2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 4224,1},{-2212,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode3CC = { {-4224,1},{ 2212,1}, { 211,1}};
      static const map<PdgId,unsigned int> & mode4   = { { 4122,1},{-2212,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & mode4CC = { {-4122,1},{ 2212,1}, { 211,1}, {-211,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign = 1,imode=-1, ibaryon,ipi=-211;
	if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode1)) {
	  sign=1;
	  imode=0;
	  ibaryon=4112;
	  ipi=211;
	}
	else if  (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode1CC)) {
	  sign=-1;
	  imode=0;
	  ibaryon=4112;
	  ipi=211;
	}
	else if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode2)) {
	  sign=1;
	  imode=1;
	  ibaryon=4222;
	}
	else if  (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode2CC)) {
	  sign=-1;
	  imode=1;
	  ibaryon=4222;
	}
	else if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode3)) {
	  sign=1;
	  imode=2;
	  ibaryon=4224;
	}
	else if  (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode3CC)) {
	  sign=-1;
	  imode=2;
	  ibaryon=4224;
	}
	else if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,4,mode4)) {
	  sign=1;
	  imode=3;
	  ibaryon=4122;
	}
	else if  (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,4,mode4CC)) {
	  sign=-1;
	  imode=3;
	  ibaryon=4122;
	}
	else
	  continue;
	const Particle & SigC = B0.decayProducts()[ix].at( sign*ibaryon)[0];
       	const Particle & pbar = B0.decayProducts()[ix].at(-sign*2212   )[0];
       	const Particle & pim  = B0.decayProducts()[ix].at( sign*ipi    )[0];
	if(imode<=2) {
	  _h[imode][0]->fill((SigC.momentum()+pim .momentum()).mass());
	  _h[imode][1]->fill((SigC.momentum()+pbar.momentum()).mass());
	  _h[imode][2]->fill((pbar.momentum()+pim .momentum()).mass());
	}
	else {
	  const Particle & pip  = B0.decayProducts()[ix].at(-sign*ipi    )[0];
	  _h[imode][0]->fill((SigC.momentum()+pip .momentum()).mass());
	  _h[imode][1]->fill((SigC.momentum()+pim .momentum()).mass());
	  _h[imode][2]->fill((pbar.momentum()+pip .momentum()).mass());
	  _h[imode][3]->fill((pbar.momentum()+pim .momentum()).mass());
	  FourMomentum ppipi = pip .momentum()+pim .momentum();
	  _h[imode][4]->fill(ppipi.mass());
	  FourMomentum plam = SigC.momentum()+pbar.momentum();
	  _h[imode][5]->fill(plam.mass());
	  _h[4][0]->fill((ppipi+pbar.momentum()).mass());
	  _h[4][1]->fill((ppipi+SigC.momentum()).mass());
	  _h[4][2]->fill((plam+pim.momentum()).mass());
	  _h[4][3]->fill((plam+pip.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix) {
	for(unsigned int iy=0;iy<6;++iy) {
	  if( (ix<3 && iy>=3) or (ix==4 && iy>=5) ) continue;
	  normalize(_h[ix][iy],1.,false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5][6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2013_I1217425);

}
