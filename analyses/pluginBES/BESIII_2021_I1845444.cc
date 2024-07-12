// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  Ds -> KS0 K-pi+pi+
  class BESIII_2021_I1845444 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1845444);


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
      for(unsigned int ix=0;ix<9;++ix)
	book(_h[ix],1,1,1+ix);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2}, {-321,1}, {310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,2}, { 321,1}, {310,1}};
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
	const Particle  & Km  = DS.decayProducts()[ix].at(-sign*321)[0];
	const Particles & pip = DS.decayProducts()[ix].at( sign*211);
	const Particle  & K0  = DS.decayProducts()[ix].at(      310)[0];
	double mK0pi[2] = {(K0.momentum()+pip[0].momentum()).mass(), 
			   (K0.momentum()+pip[1].momentum()).mass()};
	unsigned int i0(0),i1(1);
	if(mK0pi[i0]>mK0pi[i1]) swap(i0,i1);
	_h[0]->fill((K0 .momentum()+Km .momentum()).mass());
	_h[1]->fill(mK0pi[i0]);
	_h[2]->fill((Km .momentum()+pip[i1].momentum()).mass());
	_h[3]->fill(mK0pi[i1]);
	_h[4]->fill((Km .momentum()+pip[i0].momentum()).mass());
	_h[5]->fill((K0 .momentum()+Km .momentum()+pip[i0].momentum()).mass());
	_h[6]->fill((K0 .momentum()+Km .momentum()+pip[i1].momentum()).mass());
	_h[7]->fill((pip[i0].momentum()+pip[i1].momentum()).mass());
	_h[8]->fill((Km .momentum()+pip[i0].momentum()+pip[i1].momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1845444);

}
