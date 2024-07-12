// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+pi+pi-
  class BABAR_2008_I792597 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I792597);


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
      book(_h_pippim[0],1,1,3);
      book(_h_pippip   ,1,1,4);
      book(_h_pippim[1],1,1,1);
      book(_h_pippim[2],1,1,2);
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
	// kinematic variables
	double x[3] = {(pim.momentum()+pip[0].momentum()).mass2(),
		       (pim.momentum()+pip[1].momentum()).mass2(),
		       (pip[0].momentum()+pip[1].momentum()).mass2()};
	if(x[0]>x[1]) swap(x[0],x[1]);
	_dalitz->fill(x[0],x[1]);
	_dalitz->fill(x[1],x[0]);
	// fill plots
	_h_pippim[0]->fill(x[0]);
	_h_pippim[0]->fill(x[1]);
	_h_pippim[1]->fill(x[0]);
	_h_pippim[2]->fill(x[1]);
	_h_pippip   ->fill(x[2]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_pippim[ix]);
      normalize(_h_pippip);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim[3],_h_pippip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I792597);

}
