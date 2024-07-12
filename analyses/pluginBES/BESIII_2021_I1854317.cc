// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> KS0 pi+ pi0
  class BESIII_2021_I1854317 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1854317);


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
      declare(DS, "DS");
      // histograms
      book(_h_Kpi0,1,1,1);
      book(_h_Kpip,1,1,2);
      book(_h_pipi,1,1,3);
      book(_dalitz,"dalitz",50,0.3,3.5,50,0.3,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 111,1}, {310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 111,1}, {310,1}};
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
	const Particle & K0  = DS.decayProducts()[ix].at(      310)[0];
	const Particle & pi0 = DS.decayProducts()[ix].at(      111)[0];
	const Particle & pip = DS.decayProducts()[ix].at( sign*211)[0];
	double mplus = (K0 .momentum()+pip.momentum()).mass2();
	double mzero = (K0 .momentum()+pi0.momentum()).mass2();
	double mpipi = (pip.momentum()+pi0.momentum()).mass2();
	_dalitz->fill(mzero,mplus);
	_h_pipi->fill(sqrt(mpipi));
	_h_Kpip->fill(sqrt(mplus));
	_h_Kpi0->fill(sqrt(mzero));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi0);
      normalize(_h_Kpip);
      normalize(_h_pipi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi0,_h_Kpip,_h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1854317);

}
