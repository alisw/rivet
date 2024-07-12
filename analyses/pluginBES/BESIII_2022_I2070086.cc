// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> KS0 K+ pi0
  class BESIII_2022_I2070086 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2070086);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      declare(DS, "DS");
      // histograms
      book(_h_KK  ,1,1,1);
      book(_h_K0pi,1,1,2);
      book(_h_Kppi,1,1,3);
      book(_dalitz,"dalitz",50,0.3,2.3,50,0.3,2.3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{ 111,1}, {310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 111,1}, {310,1}};
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
	const Particle & Kp  = DS.decayProducts()[ix].at( sign*321)[0];
	double m0  = (K0.momentum()+pi0.momentum()).mass2();
	double mp  = (Kp.momentum()+pi0.momentum()).mass2();
	double mKK = (K0.momentum()+Kp.momentum()).mass2();
	_dalitz->fill(mp,m0);
	_h_KK->fill(sqrt(mKK));
	_h_Kppi->fill(sqrt(mp));
	_h_K0pi->fill(sqrt(m0));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_K0pi);
      normalize(_h_Kppi);
      normalize(_h_KK  );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0pi,_h_Kppi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2070086);

}
