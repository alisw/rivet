// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D -> eta' decays
  class BESIII_2018_I1693610 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1693610);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411||
						Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      DD.addStable(PID::ETA);
      DD.addStable(PID::ETAPRIME);
      declare(DD, "DD");
      // histograms
      for(unsigned int ix=0;ix<9;++ix)
	book(_h[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<3;++ix)
	book(_dalitz[ix],"dalitz_"+toString(ix+1),50,0.4,0.9,50,1.1,1.9);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay modes
      static const map<PdgId,unsigned int> & mode1   = { {-321,1},{ 211,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 321,1},{-211,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,1},{ 111,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 310,1},{ 211,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 310,1},{-211,1}, { 331,1}};
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	// D0 -> K- pi+ omega
	if( (DD.decaying()[ix].pid()== 421 && DD.modeMatches(ix,3,mode1)) ||
	    (DD.decaying()[ix].pid()==-421 && DD.modeMatches(ix,3,mode1CC))) {
	  int sign = DD.decaying()[ix].pid()/421;
	  const Particles & pip  = DD.decayProducts()[ix].at( sign*211);
	  const Particles & Km   = DD.decayProducts()[ix].at(-sign*321);
	  const Particles & etaP = DD.decayProducts()[ix].at(331);
	  double mKpi    = (pip[0].momentum()+Km[0].momentum()).mass2();
	  double mpietaP = (pip[0].momentum()+etaP[0].momentum()).mass2();
	  _dalitz[0]->fill(mKpi,mpietaP);
	  _h[0]->fill(sqrt(mKpi));
	  _h[1]->fill((Km[0].momentum()+etaP[0].momentum()).mass());
	  _h[2]->fill(sqrt(mpietaP));
	}
	// D0 -> KS0 pi0 etaP
	else if (DD.decaying()[ix].abspid()==421 && DD.modeMatches(ix,3,mode2)) {
	  const Particles & pi0  = DD.decayProducts()[ix].at(111);
	  const Particles & KS0  = DD.decayProducts()[ix].at(310);
	  const Particles & etaP = DD.decayProducts()[ix].at(331);
	  double mKpi    = (pi0[0].momentum()+KS0[0].momentum()).mass2();
	  double mpietaP = (pi0[0].momentum()+etaP[0].momentum()).mass2();
	  _dalitz[1]->fill(mKpi,mpietaP);
	  _h[3]->fill(sqrt(mKpi));
	  _h[4]->fill((KS0[0].momentum()+etaP[0].momentum()).mass());
	  _h[5]->fill(sqrt(mpietaP));
	}
	// D0 -> K- pi+ etaP
	else if( (DD.decaying()[ix].pid()== 411 && DD.modeMatches(ix,3,mode3)) ||
		 (DD.decaying()[ix].pid()==-411 && DD.modeMatches(ix,3,mode3CC))) {
	  int sign = DD.decaying()[ix].pid()/411;
	  const Particles & pip  = DD.decayProducts()[ix].at( sign*211);
	  const Particles & KS0  = DD.decayProducts()[ix].at(310);
	  const Particles & etaP = DD.decayProducts()[ix].at(331);
	  double mKpi    = (pip[0].momentum()+KS0[0].momentum()).mass2();
	  double mpietaP = (pip[0].momentum()+etaP[0].momentum()).mass2();
	  _dalitz[2]->fill(mKpi,mpietaP);
	  _h[6]->fill(sqrt(mKpi));
	  _h[7]->fill((KS0[0].momentum()+etaP[0].momentum()).mass());
	  _h[8]->fill(sqrt(mpietaP));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h[ix],1.,false);
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_dalitz[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    Histo2DPtr _dalitz[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1693610);

}
