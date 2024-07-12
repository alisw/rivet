// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  chi_c -> pbar Sigma0 K+
  class BESIII_2020_I1818254 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1818254);


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
      chi.addStable( PID::SIGMA0);
      chi.addStable(-PID::SIGMA0);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<3;++ix) {
	book(_dalitz[ix], "dalitz_"+toString(ix+1),50,2.5,7.,50,1.8,6.);
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-2212,1}, { 3212,1}, { 321,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 2212,1}, {-3212,1}, {-321,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	int sign=1;
	if(chi.modeMatches(ix,3,mode)) {
	  sign =  1;
	}
	else if(chi.modeMatches(ix,3,modeCC)) {
	  sign = -1;
	}
	else continue;
	unsigned int iloc = chi.decaying()[ix].pid()==10441 ? 0 : chi.decaying()[ix].pid()==445 ? 2 : 1;
	const Particle & Kp   = chi.decayProducts()[ix].at( sign*321)[0];
	const Particle & pbar = chi.decayProducts()[ix].at(-sign*2212)[0];
	const Particle & sig  = chi.decayProducts()[ix].at(sign*3212)[0];
	double mpK   = (Kp .momentum()+pbar.momentum()).mass2();
	double msigK = (Kp .momentum()+sig .momentum()).mass2();
	double msigp = (sig.momentum()+pbar.momentum()).mass2();
	_h[iloc][0]->fill(sqrt(mpK));
	_h[iloc][1]->fill(sqrt(msigK));
	_h[iloc][2]->fill(sqrt(msigp));
	_dalitz[iloc]->fill(msigK,mpK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_dalitz[ix]);
	for(unsigned int iy=0;iy<3;++iy) {
	  normalize(_h[ix][iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][3];
    Histo2DPtr _dalitz[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1818254);

}
