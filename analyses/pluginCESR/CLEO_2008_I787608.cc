// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @briefchi_c -> pi+pi- 2pi0 and K pi K0 pi0
  class CLEO_2008_I787608 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2008_I787608);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441 or
						Cuts::pid==20443 or
						Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::ETA);
      chi.addStable( PID::ETAPRIME);
      declare(chi, "chi");
      unsigned int iK[3]={4,2,5};
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<iK[ix];++iy) {
	  // 4 pion
	  if(iy<2) book(_h_pi[ix][iy],1+ix,1,1+iy);
	  // 2 K 2 pi
	  book(_h_K[ix][iy],4+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 211,1}, { -211,1}, {111,2} };
      static const map<PdgId,unsigned int> & mode2 = { { 321,1}, { -211,1}, {111,1} , {310,1}};
      static const map<PdgId,unsigned int> & mode3 = { {-321,1}, {  211,1}, {111,1} , {310,1}};
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	unsigned int iloc = 0;
	int sign=1;
	if   (chi.decaying()[ix].pid()==20443) iloc=1;
	else if(chi.decaying()[ix].pid()==445) iloc=2;
	if(chi.modeMatches(ix,4,mode1)) {
	  const Particle & pim  = chi.decayProducts()[ix].at(-211)[0];
	  const Particle & pip  = chi.decayProducts()[ix].at( 211)[0];
	  const Particles & pi0 = chi.decayProducts()[ix].at( 111);
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h_pi[iloc][0]->fill((pip.momentum()+pi0[ix].momentum()).mass());
	    _h_pi[iloc][1]->fill((pim.momentum()+pi0[ix].momentum()).mass());
	  }
	  continue;
	}
	else if(chi.modeMatches(ix,4,mode2)) {
	  sign=1;
	}
	else if(chi.modeMatches(ix,4,mode3)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & pim = chi.decayProducts()[ix].at(-sign*211)[0];
	const Particle & Kp  = chi.decayProducts()[ix].at( sign*321)[0];
	const Particle & pi0 = chi.decayProducts()[ix].at( 111)[0];
	const Particle & K0  = chi.decayProducts()[ix].at( 310)[0];
	_h_K[iloc][0]->fill((Kp .momentum()+pim.momentum()).mass());
	_h_K[iloc][1]->fill((pim.momentum()+pi0.momentum()).mass());
	if (_h_K[iloc][2])
	  _h_K[iloc][2]->fill((Kp.momentum()+pi0.momentum()).mass());
	if (_h_K[iloc][3])
	  _h_K[iloc][3]->fill((K0.momentum()+pim.momentum()).mass());
	if (_h_K[iloc][4])
	  _h_K[iloc][4]->fill((K0.momentum()+pi0.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      unsigned int iK[3]={4,2,5};
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<iK[ix];++iy) {
	  // 4 pion
	  if(iy<2) normalize(_h_pi[ix][iy]);
	  // 2 K 2 pi
	  normalize(_h_K[ix][iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi[3][2];
    Histo1DPtr _h_K [3][5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2008_I787608);

}
