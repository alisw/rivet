// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> phi phi eta
  class BESIII_2020_I1763897 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1763897);


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
      chi.addStable( PID::ETAPRIME);
      chi.addStable( PID::PHI);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<3;++ix) {
	book(_dalitz[ix], "dalitz_"+toString(ix+1),50,2.,7.,50,2.,7.);
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {333,2}, {221,1}};
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,3,mode)) continue;
	unsigned int iloc = chi.decaying()[ix].pid()==10441 ? 0 : chi.decaying()[ix].pid()==445 ? 2 : 1;
	const Particle  & eta = chi.decayProducts()[ix].at(221)[0];
	const Particles & phi = chi.decayProducts()[ix].at(333);
	_h[iloc][0]->fill((phi[0].momentum()+phi[1].momentum()).mass());
	double mPhiEta[2]={(phi[0].momentum()+eta.momentum()).mass2(),
			   (phi[1].momentum()+eta.momentum()).mass2()};
	_h[iloc][1]->fill(sqrt(mPhiEta[0]));
	_h[iloc][1]->fill(sqrt(mPhiEta[1]));
	_dalitz[iloc]->fill(mPhiEta[0],mPhiEta[1]);
	_dalitz[iloc]->fill(mPhiEta[1],mPhiEta[0]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_dalitz[ix]);
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h[ix][iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    Histo2DPtr _dalitz[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1763897);

}
