// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  Ds -> K+pi+pi-pi0
  class BESIII_2022_I2088218 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2088218);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      declare(DS,"DS");
      // histos
      for(unsigned int ix=0;ix<10;++ix)
	book(_h[ix],1,1,1+ix);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 321,1}, {111,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 211,1}, {-211,1}, {-321,1}, {111,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	int sign = 1;
	if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & pip = DS.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim = DS.decayProducts()[ix].at(-sign*211)[0];
	const Particle & Kp  = DS.decayProducts()[ix].at( sign*321)[0];
	const Particle & pi0 = DS.decayProducts()[ix].at(      111)[0];
	double mpipi = (pim.momentum()+pip.momentum()).mass();
	if (mpipi>.46 && mpipi<.52) continue;
	_h[0]->fill((Kp .momentum()+pim.momentum()).mass());
	_h[1]->fill((Kp .momentum()+pi0.momentum()).mass());
	_h[2]->fill(mpipi);
	_h[3]->fill((pip.momentum()+pi0.momentum()).mass());
	_h[4]->fill((pim.momentum()+pi0.momentum()).mass());
	_h[5]->fill((Kp .momentum()+pim.momentum()+pi0.momentum()).mass());
	_h[6]->fill((pip.momentum()+pim.momentum()+pi0.momentum()).mass());
	_h[7]->fill((Kp .momentum()+pip.momentum()).mass());
	_h[8]->fill((Kp .momentum()+pip.momentum()+pi0.momentum()).mass());
	_h[9]->fill((Kp .momentum()+pip.momentum()+pim.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<10;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[10];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2088218);

}
