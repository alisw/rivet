// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_cJ -> 4KS0
  class BESIII_2019_I1716627 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1716627);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443 or
						Cuts::pid==445   or
						Cuts::pid==10441);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::ETA);
      declare(chi, "chi");
      for(unsigned int iy=0;iy<2;++iy) {
	book(_h[iy],1,1,1+iy);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {310,4} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
       	if(!chi.modeMatches(ix,4,mode)) continue;
	const Particles & K0   = chi.decayProducts()[ix].at(310);
	for(unsigned int iy=0;iy<4;++iy) {
	  _h[1]->fill((chi.decaying()[ix].momentum()-K0[iy].momentum()).mass());
	  for(unsigned int iz=ix+1;iz<4;++iz)
	    _h[0]->fill((K0[ix].momentum()+K0[iz].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<2;++iy) {
	normalize(_h[iy]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1716627);

}
