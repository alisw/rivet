// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+pi+pi-
  class E791_2001_I530319 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E791_2001_I530319);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      declare(DS,"DS");
      // histos
      book(_h_pi,1,1,1);
      book(_dalitz, "dalitz",50,0.,3.5,50,0.0,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,2}, { 211,1}};
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
	const Particles & pip = DS.decayProducts()[ix].at( sign*211);
	const Particle  & pim = DS.decayProducts()[ix].at(-sign*211)[0];
	double m1 = (pim.momentum()+pip[0].momentum()).mass2();
	double m2 = (pim.momentum()+pip[1].momentum()).mass2();
	_h_pi->fill(m1);
	_h_pi->fill(m2);
	_dalitz->fill(m1,m2);
	_dalitz->fill(m2,m1);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E791_2001_I530319);

}
