// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  chi_c1 -> eta pi+pi-
  class BESIII_2017_I1490896 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1490896);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::ETA);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,ix+1);
      book(_dalitz,"dalitz",50,0.,12.,50,0., 9.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 221,1}, { 211,1}, {-211,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(chi.modeMatches(ix,3,mode1)) {
	  const Particle & eta = chi.decayProducts()[ix].at( 221)[0];
	  const Particle & pim = chi.decayProducts()[ix].at(-211)[0];
	  const Particle & pip = chi.decayProducts()[ix].at( 211)[0];
	  double m1 = (pip.momentum()+pim.momentum()).mass2();
	  double m2 = (eta.momentum()+pip.momentum()).mass2();
	  double m3 = (eta.momentum()+pim.momentum()).mass2();
	  _dalitz->fill(m2,m1);
	  _h[1]->fill(sqrt(m1));
	  _h[0]->fill(sqrt(m2));
	  _h[0]->fill(sqrt(m3));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_dalitz);
      for(unsigned int iy=0;iy<2;++iy)
	normalize(_h[iy]);
    }
    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    Histo2DPtr _dalitz;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1490896);

}
