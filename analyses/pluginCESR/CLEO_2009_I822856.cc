// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> omega pi+ pi0
  class CLEO_2009_I822856 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2009_I822856);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      DS.addStable(PID::ETA);
      DS.addStable(PID::OMEGA);
      DS.addStable(PID::ETAPRIME);
      declare(DS, "DS");
      // histograms
      book(_h_pipi,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 111,1}, { 223,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 111,1}, { 223,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	int sign = 1;
	if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & pi0 = DS.decayProducts()[ix].at(      111);
	const Particles & pip = DS.decayProducts()[ix].at( sign*211);
	_h_pipi->fill((pi0[0].momentum()+pip[0].momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2009_I822856);

}
