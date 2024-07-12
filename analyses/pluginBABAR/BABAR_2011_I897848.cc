// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> K+ pi- pi0
  class BABAR_2011_I897848 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2011_I897848);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::PI0);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 111,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 111,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & Kp  = B0.decayProducts()[ix].at( sign*321)[0];
	const Particle & pim = B0.decayProducts()[ix].at(-sign*211)[0];
	const Particle & pi0 = B0.decayProducts()[ix].at(      111)[0];
	double mpipi = (pim.momentum()+pi0.momentum()).mass();
	_h[0][0]->fill(mpipi);
	if(mpipi<1.8) _h[0][1]->fill(mpipi);
	double mKpi  = (Kp.momentum()+pi0.momentum()).mass();
	_h[1][0]->fill(mKpi);
	if(mKpi<1.8) _h[1][1]->fill(mKpi);
	mKpi  = (Kp.momentum()+pim.momentum()).mass();
	_h[2][0]->fill(mKpi);
	if(mKpi<1.8) _h[2][1]->fill(mKpi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h[ix][iy],1.,false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2011_I897848);

}
