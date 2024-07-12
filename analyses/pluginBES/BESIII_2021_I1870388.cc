// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> n K0S Lambdabar +cc
  class BESIII_2021_I1870388 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1870388);


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
      chi.addStable( PID::LAMBDA);
      chi.addStable(-PID::LAMBDA);
      declare(chi, "chi");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_dalitz[ix], "dalitz_"+toString(ix+1),50,2.,7.,50,2.,6.2);
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h[ix][iy],ix+1,1,iy+1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2112,1}, {-3122,1}, { 310,1} };
      static const map<PdgId,unsigned int> & modeCC = { {-2112,1}, { 3122,1}, { 310,1} };
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
	const Particle & KS0  = chi.decayProducts()[ix].at( 310)[0];
	const Particle & nn   = chi.decayProducts()[ix].at( sign*2112)[0];
	const Particle & lbar = chi.decayProducts()[ix].at(-sign*3122)[0];
	double mnLam = (nn  .momentum()+lbar.momentum()).mass2();
	double mnK   = (nn  .momentum()+ KS0.momentum()).mass2();
	double mKLam = (lbar.momentum()+ KS0.momentum()).mass2();
	_h[iloc][0]->fill(sqrt(mnLam));
	_h[iloc][1]->fill(sqrt(mnK  ));
	_h[iloc][2]->fill(sqrt(mKLam));
	_dalitz[iloc]->fill(mKLam,mnK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_dalitz[ix]);
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h[ix][iy]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][3];
    Histo2DPtr _dalitz[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1870388);

}
