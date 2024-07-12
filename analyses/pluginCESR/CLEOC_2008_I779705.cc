// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> eta pi+pi- decays
  class CLEOC_2008_I779705 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2008_I779705);


    /// @name Analysis methods
    ///@{

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
      book(_h_eta_pi,1,1,1);
      book(_h_pi_pi ,1,1,2);
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{-211,1}, {221,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( !D0.modeMatches(ix,3,mode)  ) continue;
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & eta = D0.decayProducts()[ix].at(221);
	const Particles & pip = D0.decayProducts()[ix].at( sign*211);
	const Particles & pim = D0.decayProducts()[ix].at(-sign*211);
	_h_eta_pi->fill((pip[0].momentum()+eta[0].momentum()).mass());
	_h_pi_pi ->fill((pip[0].momentum()+pim[0].momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_eta_pi);
      normalize(_h_pi_pi );
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_eta_pi,_h_pi_pi;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I779705);

}
