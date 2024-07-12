// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0, KL0 pi+ pi-
  class BESIII_2022_I2615968 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2615968);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // Histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // // define the decay mode
      static const map<PdgId,unsigned int> & mode1 = { { 310,1}, { 211,1},{-211,1}};
      static const map<PdgId,unsigned int> & mode2 = { { 130,1}, { 211,1},{-211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
       	int sign = D0.decaying()[ix].pid()/421;
	int imode = 0;
       	// KS0 pi+pi-
       	if      (D0.modeMatches(ix,3,mode1) ) imode=0;
	else if (D0.modeMatches(ix,3,mode2) ) imode=1;
	else continue;
       	const Particle & pip= D0.decayProducts()[ix].at( sign*211)[0];
       	const Particle & pim= D0.decayProducts()[ix].at(-sign*211)[0];
       	const Particle & K0 = D0.decayProducts()[ix].at( imode==0 ? 310 : 130)[0];
      	double mminus = (pim.momentum()+K0.momentum() ).mass2();
      	double mplus  = (pip.momentum()+K0.momentum() ).mass2();
      	double mpipi  = (pip.momentum()+pim.momentum()).mass2();
       	_h[imode][0]->fill(mplus );
       	_h[imode][1]->fill(mminus);
       	_h[imode][2]->fill(mpipi );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2615968);

}
