// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> pbar K+ Sigma^$ and chi_cJ -> pbar K+Lambda 
  class BESIII_2013_I1203840 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1203840);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443 or
						Cuts::pid==445   or
						Cuts::pid==10441 or
						Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::SIGMA0);
      chi.addStable( PID::SIGMA0);
      chi.addStable( PID::LAMBDA);
      chi.addStable(-PID::LAMBDA);
      declare(chi, "chi");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_psi[ix],1,1,ix+1);
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h_chi[iy][ix],2,iy+1,ix+1);
      }
      book(_h_pLam,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { {-2212,1}, { 3212,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 2212,1}, {-3212,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode2   = { {-2212,1}, { 3122,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode2CC = { { 2212,1}, {-3122,1}, {-321,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
      	int sign=1;
	if(chi.decaying()[ix].pid()==100443) {
	  if(chi.modeMatches(ix,3,mode1)) {
	    sign =  1;
	  }
	  else if(chi.modeMatches(ix,3,mode1CC)) {
	    sign = -1;
	  }
	  else continue;
	  const Particle & pbar  = chi.decayProducts()[ix].at(-sign*2212)[0];
	  const Particle & Kp    = chi.decayProducts()[ix].at( sign*321 )[0];
	  const Particle & sigma = chi.decayProducts()[ix].at( sign*3212)[0];
	  _h_psi[0]->fill((pbar.momentum()+sigma.momentum()).mass());
	  _h_psi[1]->fill((Kp  .momentum()+sigma.momentum()).mass());
	}
	else {
	  if(chi.modeMatches(ix,3,mode2)) {
	    sign =  1;
	  }
	  else if(chi.modeMatches(ix,3,mode2CC)) {
	    sign = -1;
	  }
	  else continue;
	  unsigned int iloc = chi.decaying()[ix].pid()==10441 ? 0 : chi.decaying()[ix].pid()==445 ? 2 : 1;
	  const Particle & pbar = chi.decayProducts()[ix].at(-sign*2212)[0];
	  const Particle & Kp   = chi.decayProducts()[ix].at( sign*321 )[0];
	  const Particle & lam  = chi.decayProducts()[ix].at( sign*3122)[0];
	  _h_chi[iloc][0]->fill((pbar.momentum()+Kp .momentum()).mass());
	  _h_chi[iloc][1]->fill((Kp  .momentum()+lam.momentum()).mass());
	  if(iloc==0) _h_pLam->fill((pbar.momentum()+lam.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_psi[ix],1.,false);
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h_chi[iy][ix],1.,false);
      }
      normalize(_h_pLam,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_psi[2],_h_chi[3][2],_h_pLam;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1203840);

}
