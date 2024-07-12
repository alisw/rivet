// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c0 -> pi+ pi- K+ K-
  class BESII_2005_I690784 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2005_I690784);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::ETA);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,ix+1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 211,1}, {-211,1}, {321,1}, {-321,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,4,mode1)) continue;
	const Particle & pim = chi.decayProducts()[ix].at(-211)[0];
	const Particle & pip = chi.decayProducts()[ix].at( 211)[0];
	const Particle & Km  = chi.decayProducts()[ix].at(-321)[0];
	const Particle & Kp  = chi.decayProducts()[ix].at( 321)[0];
	_h[0]->fill((Km .momentum() +Kp.momentum()).mass());
	_h[1]->fill((pim.momentum()+pip.momentum()).mass());
	_h[2]->fill((Km .momentum()+pip.momentum()).mass());
	_h[2]->fill((Kp .momentum()+pim.momentum()).mass());
	_h[3]->fill((Km .momentum()+pip.momentum()+pim.momentum()).mass());
	_h[3]->fill((Kp .momentum()+pim.momentum()+pip.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2005_I690784);

}
