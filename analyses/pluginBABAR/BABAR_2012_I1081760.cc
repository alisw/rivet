// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> pi+ pi- K*0
  class BABAR_2012_I1081760 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1081760);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 313);
      B0.addStable(-313);
      declare(B0, "B0");
      // histogram
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 313,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 211,1}, {-211,1}, {-313,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode))        sign= 1;
      	else if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) sign=-1;
	else continue;
	const Particle & pip = B0.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim = B0.decayProducts()[ix].at(-sign*211)[0];
	_h->fill((pip.momentum()+pim.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1081760);

}
