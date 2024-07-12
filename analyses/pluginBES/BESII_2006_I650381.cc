// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> p nbar pi+ + c.c
  class BESII_2006_I650381 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I650381);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix],1,1,1+ix);
	book(_dalitz[ix], "dalitz_"+toString(ix+1),50,1.,5.,50,1.,5.);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2112,1}, {-211,1} };
      static const map<PdgId,unsigned int> & modeCC = { {-2212,1}, { 2112,1}, { 211,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(psi.modeMatches(ix,3,mode)) {
	  const Particle & pim  = psi.decayProducts()[ix].at(-211)[0];
	  const Particle & pp   = psi.decayProducts()[ix].at( 2212)[0];
	  const Particle & nbar = psi.decayProducts()[ix].at(-2112)[0];
	  double mminus = (pp  .momentum()+pim.momentum()).mass2();
	  double mplus  = (nbar.momentum()+pim.momentum()).mass2();
	  _h[0]->fill(sqrt(mminus));
	  _dalitz[0]->fill(mminus,mplus);
	}
	else if (psi.modeMatches(ix,3,modeCC)) {
	  const Particle & pip  = psi.decayProducts()[ix].at( 211)[0];
	  const Particle & nn   = psi.decayProducts()[ix].at( 2112)[0];
	  const Particle & pbar = psi.decayProducts()[ix].at(-2212)[0];
	  double mminus = (pbar.momentum()+pip.momentum()).mass2();
	  double mplus  = (nn  .momentum()+pip.momentum()).mass2();
	  _h[1]->fill(sqrt(mminus));
	  _dalitz[1]->fill(mminus,mplus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h[ix]);
	normalize(_dalitz[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I650381);

}
