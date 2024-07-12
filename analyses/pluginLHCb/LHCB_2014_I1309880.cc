// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_c -> J/psi p pbar pi+
  class LHCB_2014_I1309880 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2014_I1309880);


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
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 443,1}, { 2212,1}, { 2212,1}, { 211,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 443,1}, { 2212,1}, { 2212,1}, {-211,1} };
      DecayedParticles BC = apply<DecayedParticles>(event, "BC");
      // loop over particles
      for(unsigned int ix=0;ix<BC.decaying().size();++ix) {
	int sign = BC.decaying()[ix].pid()/BC.decaying()[ix].abspid();
	if ((sign== 1 && BC.modeMatches(ix,4,mode  )) ||
	    (sign==-1 && BC.modeMatches(ix,4,modeCC))) {
	  const Particle & prot = BC.decayProducts()[ix].at( sign*2212)[0];
	  const Particle & pbar = BC.decayProducts()[ix].at(-sign*2212)[0];
	  const Particle & pip  = BC.decayProducts()[ix].at( sign*211 )[0];
	  _h[0]->fill((prot.momentum()+pbar.momentum()).mass());
	  _h[1]->fill((prot.momentum()+pip .momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2014_I1309880);

}
