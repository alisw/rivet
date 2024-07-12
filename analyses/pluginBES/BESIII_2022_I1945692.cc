// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> KS0 KS0 pi+
  class BESIII_2022_I1945692 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I1945692);


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
      // histograms
      book(_h_KK ,1,1,1);
      book(_h_Kpi,1,1,2);
      book(_dalitz,"dalitz",50,0.3,2.3,50,0.3,2.3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {310,2}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1}, {310,2}};
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
	const Particles & K0  = DS.decayProducts()[ix].at(      310);
	const Particle  & pip = DS.decayProducts()[ix].at( sign*211)[0];
	double m1  = (K0[0].momentum()+pip.momentum()).mass2();
	double m2  = (K0[1].momentum()+pip.momentum()).mass2();
	double mKK = (K0[0].momentum()+K0[1].momentum()).mass2();
	_dalitz->fill(m1,m2);
	_dalitz->fill(m2,m1);
	_h_KK->fill(sqrt(mKK));
	_h_Kpi->fill(sqrt(m1));
	_h_Kpi->fill(sqrt(m2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi);
      normalize(_h_KK );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I1945692);

}
