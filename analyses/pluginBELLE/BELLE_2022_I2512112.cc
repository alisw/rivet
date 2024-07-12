// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D pi (pi) semileptonic
  class BELLE_2022_I2512112 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2512112);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 or
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(411); BB.addStable(-411);
      BB.addStable(421); BB.addStable(-421);
      BB.addStable(413); BB.addStable(-413);
      BB.addStable(423); BB.addStable(-423);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> mode1[4] = {{ {-421,1}, {-211,1}, {-11,1}, { 12,1}},
						       { { 421,1}, { 211,1}, { 11,1}, {-12,1}},
						       { {-421,1}, {-211,1}, {-13,1}, { 14,1}},
						       { { 421,1}, { 211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode2[4] = {{ {-411,1}, { 211,1}, {-11,1}, { 12,1}},
						       { { 411,1}, {-211,1}, { 11,1}, {-12,1}},
						       { {-411,1}, { 211,1}, {-13,1}, { 14,1}},
						       { { 411,1}, {-211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode3[4] = {{ {-423,1}, {-211,1}, {-11,1}, { 12,1}},
						       { { 423,1}, { 211,1}, { 11,1}, {-12,1}},
						       { {-423,1}, {-211,1}, {-13,1}, { 14,1}},
						       { { 423,1}, { 211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode4[4] = {{ {-413,1}, { 211,1}, {-11,1}, { 12,1}},
						       { { 413,1}, {-211,1}, { 11,1}, {-12,1}},
						       { {-413,1}, { 211,1}, {-13,1}, { 14,1}},
						       { { 413,1}, {-211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode5[4] = {{ {-411,1}, {-211,1}, { 211,1}, {-11,1}, { 12,1}},
						       { { 411,1}, {-211,1}, { 211,1}, { 11,1}, {-12,1}},
						       { {-411,1}, {-211,1}, { 211,1}, {-13,1}, { 14,1}},
						       { { 411,1}, {-211,1}, { 211,1}, { 13,1}, {-14,1}}};
      static const map<PdgId,unsigned int> mode6[4] = {{ {-421,1}, {-211,1}, { 211,1}, {-11,1}, { 12,1}},
						       { { 421,1}, {-211,1}, { 211,1}, { 11,1}, {-12,1}},
						       { {-421,1}, {-211,1}, { 211,1}, {-13,1}, { 14,1}},
						       { { 421,1}, {-211,1}, { 211,1}, { 13,1}, {-14,1}}};
      // loop over B mesons
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	for(unsigned int il=0;il<4;++il) {
	  int iD,iloc1,iloc2;
	  if ( BB.modeMatches(ix,4,mode1[il]) ) {
	    iloc1=0;
	    iloc2=0;
	    iD=-421;
	  }
	  else if ( BB.modeMatches(ix,4,mode2[il]) ) {
	    iloc1=0;
	    iloc2=1;
	    iD=-411;
	  }
	  else if ( BB.modeMatches(ix,4,mode3[il]) ) {
	    iloc1=1;
	    iloc2=0;
	    iD=-423;
	  }
	  else if ( BB.modeMatches(ix,4,mode4[il]) ) {
	    iloc1=1;
	    iloc2=1;
	    iD=-413;
	  }
	  else if ( BB.modeMatches(ix,5,mode5[il]) ) {
	    iloc1=2;
	    iloc2=0;
	    iD=-411;
	  }
	  else if ( BB.modeMatches(ix,5,mode6[il]) ) {
	    iloc1=2;
	    iloc2=1;
	    iD=-421;
	  }
	  else continue;
	  int sign = il%2==0 ? 1 : -1;
	  int ipi = -sign*211;
	  if(iloc2==1) ipi *=-1;
	  const Particle & pi1= BB.decayProducts()[ix].at( ipi)[0];
	  const Particle & DD = BB.decayProducts()[ix].at( iD*sign )[0];
	  FourMomentum pHad = pi1.momentum()+DD.momentum();
	  if(iloc1==2) {
	    const Particle & pi2= BB.decayProducts()[ix].at(-ipi)[0];
	    pHad += pi2.momentum();
	  }
	  double mass = pHad.mass();
	  if(iloc1==1) mass -=DD.mass();
	  _h[iloc1][iloc2]->fill(mass);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2512112);

}
