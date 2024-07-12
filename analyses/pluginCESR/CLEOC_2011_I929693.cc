// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c1 -> eta(') pi+-pi
  class CLEOC_2011_I929693 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2011_I929693);


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
      chi.addStable( PID::ETAPRIME);
      declare(chi, "chi");
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],ix+1,1,iy+1);
      book(_dalitz[0],"dalitz_1",50,0.,12.,50,0., 9.);
      book(_dalitz[1],"dalitz_2",50,0.,12.,50,0., 7.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 221,1}, { 211,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2 = { { 331,1}, { 211,1}, {-211,1} };
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
	  _dalitz[0]->fill(m2,m1);
	  _h[0][1]->fill(sqrt(m1));
	  _h[0][0]->fill(sqrt(m2));
	  _h[0][0]->fill(sqrt(m3));
	}
	else if(chi.modeMatches(ix,3,mode2)) {
	  const Particle & eta = chi.decayProducts()[ix].at( 331)[0];
	  const Particle & pim = chi.decayProducts()[ix].at(-211)[0];
	  const Particle & pip = chi.decayProducts()[ix].at( 211)[0];
	  double m1 = (pip.momentum()+pim.momentum()).mass2();
	  double m2 = (eta.momentum()+pip.momentum()).mass2();
	  double m3 = (eta.momentum()+pim.momentum()).mass2();
	  _dalitz[1]->fill(m2,m1);
	  _h[1][1]->fill(sqrt(m1));
	  _h[1][0]->fill(sqrt(m2));
	  _h[1][0]->fill(sqrt(m3));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_dalitz[ix]);
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h[ix][iy]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2011_I929693);

}
