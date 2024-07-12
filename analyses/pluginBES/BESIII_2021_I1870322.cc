// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D_s+ -> pi+ pi+ pi- eta
  class BESIII_2021_I1870322 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1870322);


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
      DS.addStable(PID::ETA);
      declare(DS,"DS");
      // histograms
      book(_h_pippip   ,1,1,1);
      book(_h_pippim   ,1,1,2);
      book(_h_pipeta   ,1,1,3);
      book(_h_pimeta   ,1,1,4);
      book(_h_3pi      ,1,1,5);
      book(_h_pippimeta,1,1,6);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2}, {-211,1}, {221,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 211,1}, {-211,2}, {221,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	int sign = 1;
	if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & pip  = DS.decayProducts()[ix].at( sign*211);
	const Particle  & pim  = DS.decayProducts()[ix].at(-sign*211)[0];
	const Particle  & eta  = DS.decayProducts()[ix].at(      221)[0];
	double metapipi[2] = {(eta.momentum()+pim.momentum()+pip[0].momentum()).mass(),
			      (eta.momentum()+pim.momentum()+pip[1].momentum()).mass()};
	if(metapipi[0]<GeV || metapipi[1]<GeV) continue;
	_h_pippip->fill((pip[0].momentum()+pip[1].momentum()).mass());
	_h_pippim->fill((pim.momentum()+pip[0].momentum()).mass());
	_h_pippim->fill((pim.momentum()+pip[1].momentum()).mass());
	_h_pipeta->fill((eta.momentum()+pip[0].momentum()).mass());
	_h_pipeta->fill((eta.momentum()+pip[1].momentum()).mass());
	_h_pimeta->fill((eta.momentum()+pim.momentum()).mass());
	_h_3pi   ->fill((pim.momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	_h_pippimeta->fill(metapipi[0]);
	_h_pippimeta->fill(metapipi[1]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pippip   );
      normalize(_h_pippim   );
      normalize(_h_pipeta   );
      normalize(_h_pimeta   );
      normalize(_h_3pi      );
      normalize(_h_pippimeta);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippip,_h_pippim,_h_pipeta,_h_pimeta,_h_3pi,_h_pippimeta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1870322);

}
