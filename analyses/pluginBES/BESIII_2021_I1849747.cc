// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Ds -> K+K-pi+pi0
  class BESIII_2021_I1849747 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1849747);


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
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, {-321,1}, { 211,1}, {111,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1}, {-321,1}, {-211,1}, {111,1}};
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
	const Particle & Kp  = DS.decayProducts()[ix].at( sign*321)[0];
	const Particle & Km  = DS.decayProducts()[ix].at(-sign*321)[0];
	const Particle & pip = DS.decayProducts()[ix].at( sign*211)[0];
	const Particle & pi0 = DS.decayProducts()[ix].at(      111)[0];
	double mKK = (Kp.momentum()+Km.momentum()).mass();
	_h[0]->fill(mKK);
	_h[1]->fill(mKK);
	_h[2]->fill((Kp .momentum()+pi0.momentum()).mass());
	_h[3]->fill((Km .momentum()+pi0.momentum()).mass());
	_h[4]->fill((pip.momentum()+pi0.momentum()).mass());
	_h[5]->fill((Km .momentum()+pip.momentum()).mass());
	_h[6]->fill((Km .momentum()+pip.momentum()+pi0.momentum()).mass());
	_h[7]->fill((Km .momentum()+ Kp.momentum()+pip.momentum()).mass());
	_h[8]->fill((Kp .momentum()+pip.momentum()+pi0.momentum()).mass());
	_h[9]->fill((Km .momentum()+ Kp.momentum()+pi0.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1849747);

}
