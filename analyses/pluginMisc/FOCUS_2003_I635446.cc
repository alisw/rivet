// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+/D+ -> pi+pi+pi-
  class FOCUS_2003_I635446 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2003_I635446);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411 or
						Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");
      // histos
      book(_h_pippim[0],1,1,1);
      book(_h_pippim[1],1,1,2);
      book(_h_pippim[2],2,1,1);
      book(_h_pippim[3],2,1,2);
      book(_dalitz[0], "dalitz1",50,0.,2.0,50,0.0,3.5);
      book(_dalitz[1], "dalitz2",50,0.,1.8,50,0.0,3.1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2},{-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,2},{ 211,1}};
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	int sign = 1;
	if (DD.decaying()[ix].pid()>0 && DD.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DD.decaying()[ix].pid()<0 && DD.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	const Particle  & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	double m1 = (pim.momentum()+pip[0].momentum()).mass2();
	double m2 = (pim.momentum()+pip[1].momentum()).mass2();
	if(m1>m2) swap(m1,m2);
	if(DD.decaying()[ix].abspid()==431) {
	  _dalitz[0]->fill(m1,m2);
	  _h_pippim[0]->fill(m1);
	  _h_pippim[1]->fill(m2);
	}
	else {
	  _dalitz[1]->fill(m1,m2);
	  _h_pippim[2]->fill(m1);
	  _h_pippim[3]->fill(m2);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_dalitz[ix]);
	normalize(_h_pippim[ix  ]);
	normalize(_h_pippim[ix+2]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim[4];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2003_I635446);

}
