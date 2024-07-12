// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> eta' K+K-
  class BESIII_2014_I1280710 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2014_I1280710);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443 or Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      chi.addStable( PID::ETA);
      chi.addStable( PID::ETAPRIME);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<4;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
      for(unsigned int ix=0;ix<2;++ix)
	book(_dalitz[ix],"dalitz_"+toString(ix+1),50,2.,10.,50,2., 10.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 331,1}, { 321,1}, {-321,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,3,mode1)) continue;
	unsigned int iloc = chi.decaying()[ix].pid()==445 ? 1 : 0;
	const Particle & eta = chi.decayProducts()[ix].at( 331)[0];
	const Particle & Km  = chi.decayProducts()[ix].at(-321)[0];
	const Particle & Kp  = chi.decayProducts()[ix].at( 321)[0];
	double m1 = (Kp .momentum()+Km.momentum()).mass2();
	double m2 = (eta.momentum()+Kp.momentum()).mass2();
	double m3 = (eta.momentum()+Km.momentum()).mass2();
	_dalitz[iloc]->fill(m2,m3);
	for(unsigned int ix=0;ix<2;++ix) {
	  _h[2*iloc+ix][0]->fill(sqrt(m1));
	  _h[2*iloc+ix][1]->fill(sqrt(m2));
	  _h[2*iloc+ix][1]->fill(sqrt(m3));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<2;++iy) {
	normalize(_dalitz[iy]);
	for(unsigned int ix=0;ix<4;++ix)
	  normalize(_h[ix][iy]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2014_I1280710);

}
