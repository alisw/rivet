// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> Lambda_c+ Lambdabar K-
  class BABAR_2011_I924163 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2011_I924163);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 3122);
      B0.addStable(-3122);
      B0.addStable( 4122);
      B0.addStable(-4122);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4122,1},{-3122,1}, {-321,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-4122,1},{ 3122,1}, { 321,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & LamC = B0.decayProducts()[ix].at( sign*4122)[0];
	const Particle & Lam  = B0.decayProducts()[ix].at(-sign*3122)[0];
	const Particle & Km   = B0.decayProducts()[ix].at(-sign*321 )[0];
	_h[0]->fill((LamC.momentum()+Lam.momentum()).mass());
	_h[1]->fill((Lam .momentum()+Km .momentum()).mass());
	_h[2]->fill((LamC.momentum()+Km .momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2011_I924163);

}
