// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+ pi0 pi0
  class BESIII_2021_I1929365 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1929365);


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
      book(_h_pi0pi0,1,1,1);
      book(_h_pippi0,1,1,2);
      book(_dalitz, "dalitz",50,0.,3.5,50,0.,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 111,2}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 111,2}};
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
	const Particle  & pip = DS.decayProducts()[ix].at( sign*211)[0];
	double mp1  = (pip   .momentum() +pi0[0].momentum()).mass2();
	double mp2  = (pip   .momentum() +pi0[1].momentum()).mass2();
	double m0   = (pi0[0].momentum() +pi0[1].momentum()).mass();
	_dalitz ->fill(mp1,mp2);
	_dalitz ->fill(mp2,mp1);
	// K_S0 veto
	if(m0>0.458 && m0<0.520) continue;
	_h_pippi0->fill(sqrt(mp1));
	_h_pippi0->fill(sqrt(mp2));
	_h_pi0pi0->fill(m0 );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pi0pi0);
      normalize(_h_pippi0);
      normalize(_dalitz  );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi0pi0,_h_pippi0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1929365);

}
