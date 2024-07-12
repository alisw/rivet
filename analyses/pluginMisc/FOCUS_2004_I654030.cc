// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+/Ds+ -> K+ pi+ pi-
  class FOCUS_2004_I654030 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2004_I654030);


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
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_Kpim[ix],1+ix,1,1);
	book(_h_pipi[ix],1+ix,1,2);
	book(_dalitz[ix],"dalitz_"+to_str(ix+1),50,0.3,3.5,50,0.,2.3);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, { 321,1}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1}, {-321,1}, { 211,1}};
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
	const Particle & Kp  = DD.decayProducts()[ix].at( sign*321)[0];
	const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	double mm    = (Kp .momentum()+pim.momentum()).mass2();
	double mpipi = (pip.momentum()+pim.momentum()).mass2();
	unsigned int iloc = DD.decaying()[ix].abspid()==411 ? 0 : 1;
	_dalitz[iloc]->fill(mm,mpipi);
	_h_pipi[iloc]->fill(mpipi);
	_h_Kpim[iloc]->fill(mm);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_Kpim[ix]);
	normalize(_h_pipi[ix]);
	normalize(_dalitz[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim[2],_h_pipi[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2004_I654030);

}
