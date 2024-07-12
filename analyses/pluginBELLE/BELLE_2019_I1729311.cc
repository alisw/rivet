// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> p pbar pi0
  class BELLE_2019_I1729311 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1729311);


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
      book(_h[0],1,1,1);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix+1],2,1,1+ix);
      book(_dalitz, "dalitz",50,1.,20.,50,1.,20.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 111,1} };
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	if(!B0.modeMatches(ix,3,mode)) continue;
	int sign = B0.decaying()[ix].pid()/511;
     	const Particle & pi0  = B0.decayProducts()[ix].at( 111)[0];
     	const Particle & pp   = B0.decayProducts()[ix].at( sign*2212)[0];
     	const Particle & pbar = B0.decayProducts()[ix].at(-sign*2212)[0];
       	double mminus = (pbar.momentum()+pi0.momentum()).mass2();
       	double mplus  = (pp  .momentum()+pi0.momentum()).mass2();
      	_h[0]->fill((pp.momentum()+pbar.momentum()).mass());
      	_h[1]->fill(sqrt(mplus ));
      	_h[2]->fill(sqrt(mminus));
       	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_dalitz);
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1729311);

}
