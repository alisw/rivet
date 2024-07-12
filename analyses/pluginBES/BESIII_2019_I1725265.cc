// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K- pi+ pi0 pi0
  class BESIII_2019_I1725265 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1725265);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { {-321,1},{ 211,1}, { 111,2}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-211,1}, { 111,2}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.decaying()[ix].pid()>0 && D0.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (D0.decaying()[ix].pid()<0 && D0.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & pi0= D0.decayProducts()[ix].at(111);
	double mpipi = (pi0[0].momentum()+pi0[1].momentum()).mass();
	if (mpipi>.458 && mpipi<.52) continue;
	const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	_h[0]->fill((Km [0].momentum()+pip[0].momentum()).mass2());
	_h[1]->fill((Km [0].momentum()+pi0[0].momentum()).mass2());
	_h[1]->fill((Km [0].momentum()+pi0[1].momentum()).mass2());
	_h[2]->fill((pip[0].momentum()+pi0[0].momentum()).mass2());
	_h[2]->fill((pip[0].momentum()+pi0[1].momentum()).mass2());
	_h[3]->fill(sqr(mpipi));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1725265);

}
