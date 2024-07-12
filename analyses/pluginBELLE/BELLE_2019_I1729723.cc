// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> KS0 K-+ pi+-
  class BELLE_2019_I1729723 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1729723);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::PI0);
      B0.addStable(PID::K0S);
      B0.addStable(PID::ETA);
      B0.addStable(PID::ETAPRIME);
      declare(B0, "B0");
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
      book(_n,"TMP/nB0");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 310,1}, {-321,1}, { 211,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 310,1}, { 321,1}, {-211,1} };
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	_n->fill();
	unsigned int iloc(0);
	int sign=1;
      	if(B0.modeMatches(ix,3,mode)) {
	  sign=1;
	  iloc = B0.decaying()[ix].pid()>0 ? 1 : 0;
	}
	else if(B0.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	  iloc = B0.decaying()[ix].pid()>0 ? 0 : 1;
	}
	else
	  continue;
     	const Particle & pip  = B0.decayProducts()[ix].at( sign*211)[0];
     	const Particle & Km   = B0.decayProducts()[ix].at(-sign*321)[0];
     	const Particle & K0   = B0.decayProducts()[ix].at(      310)[0];
	_h[0][iloc]->fill((Km.momentum()+pip.momentum()).mass());
	_h[1][iloc]->fill((K0.momentum()+pip.momentum()).mass());
	_h[2][iloc]->fill((Km.momentum()+K0 .momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h[ix][iy], 1e7/ *_n);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    CounterPtr _n;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1729723);

}
