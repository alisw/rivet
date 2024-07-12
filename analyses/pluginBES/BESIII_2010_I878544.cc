// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_cJ -> 4 pi0
  class BESIII_2010_I878544 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2010_I878544);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441 ||
						Cuts::pid==445   ||
						Cuts::pid==20443);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      declare(chi, "chi");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 111,4} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,4,mode)) continue;
	const Particles & pi0 = chi.decayProducts()[ix].at(111);
	unsigned int imode = 0;
	if(chi.decaying()[ix].pid()==20443)    imode=1;
	else if(chi.decaying()[ix].pid()==445) imode=2;
	for(unsigned int ix=0;ix<pi0.size();++ix) {
	  for(unsigned int iy=ix+1;iy<pi0.size();++iy) {
	    _h[imode]->fill((pi0[ix].momentum()+pi0[iy].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2010_I878544);

}
