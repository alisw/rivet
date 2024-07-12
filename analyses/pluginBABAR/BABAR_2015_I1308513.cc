// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> J/psi phi K
  class BABAR_2015_I1308513 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1308513);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable(333);
      BB.addStable(443);
      declare(BB, "BB");
      // histo
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 443,1},{ 333,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 443,1},{ 333,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 443,1},{ 333,1}, { 130,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 443,1},{ 333,1}, { 310,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	if (BB.modeMatches(ix,3,mode1) ||
	    BB.modeMatches(ix,3,mode1CC) ||
	    BB.modeMatches(ix,3,mode2) ||
	    BB.modeMatches(ix,3,mode2CC)) {
	  const Particle & phi = BB.decayProducts()[ix].at(333)[0];
	  const Particle & psi = BB.decayProducts()[ix].at(443)[0];
	  _h->fill((phi.momentum()+psi.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1308513);

}
