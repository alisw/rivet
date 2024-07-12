// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+ pi0 eta'
  class BESIII_2022_I2030993 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2030993);


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
      DS.addStable(PID::ETAPRIME);
      declare(DS,"DS");
      // histograms
      book(_h_etapip,1,1,1);
      book(_h_etapi0,1,1,2);
      book(_h_pipi , 1,1,3);
      book(_dalitz, "dalitz",50,0.,1.1,50,1,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, { 111,1}, { 331,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1}, { 111,1}, { 331,1}};
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
	const Particle  & pip = DS.decayProducts()[ix].at( sign*211)[0];
	const Particle  & pi0 = DS.decayProducts()[ix].at(      111)[0];
	const Particle  & eta = DS.decayProducts()[ix].at(      331)[0];
	double mplus = (eta.momentum()+pip.momentum()).mass2();
	double mzero = (eta.momentum()+pi0.momentum()).mass2();
	double mpipi = (pip.momentum()+pi0.momentum()).mass2();
	_dalitz  ->fill(mpipi,mplus);
	_h_pipi  ->fill(sqrt(mpipi));
	_h_etapip->fill(sqrt(mplus));
	_h_etapi0->fill(sqrt(mzero));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi  );
      normalize(_h_etapip);
      normalize(_h_etapi0);
      normalize(_dalitz  );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi, _h_etapip,_h_etapi0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2030993);

}
