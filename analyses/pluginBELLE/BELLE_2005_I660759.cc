// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 -> p K- K- pi+
  class BELLE_2005_I660759 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I660759);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4132);
      declare(ufs, "UFS");
      DecayedParticles XIC0(ufs);
      XIC0.addStable(PID::PI0);
      XIC0.addStable(PID::K0S);
      XIC0.addStable(PID::ETA);
      XIC0.addStable(PID::ETAPRIME);
      declare(XIC0, "XIC0");
      // histograms
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::KMINUS,2}, { PID::PROTON,1}, { PID::PIPLUS ,1}};
      static const map<PdgId,unsigned int> & modeCC = { { PID::KPLUS ,2}, {-PID::PROTON,1}, { PID::PIMINUS,1}};
      DecayedParticles XIC0 = apply<DecayedParticles>(event, "XIC0");
      // loop over particles
      for(unsigned int ix=0;ix<XIC0.decaying().size();++ix) {
	int sign = 1;
	if (XIC0.decaying()[ix].pid()>0 && XIC0.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (XIC0.decaying()[ix].pid()<0 && XIC0.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & pip = XIC0.decayProducts()[ix].at( sign*PID::PIPLUS)[0];
	const Particles & Km = XIC0.decayProducts()[ix].at( sign*PID::KMINUS);
	for(unsigned int ix=0;ix<2;++ix)
	  _h->fill((pip.momentum()+Km[ix].momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I660759);

}
