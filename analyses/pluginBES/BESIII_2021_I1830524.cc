// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D_s -> K+K-pi+
  class BESIII_2021_I1830524 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1830524);


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
      book(_h_KK[0],1,1,1);
      book(_h_KK[1],1,1,2);
      book(_h_Kmpi ,1,1,3);
      book(_h_Kppi ,1,1,4);
      book(_dalitz, "dalitz",50,0.3,3.5,50,0.07,2.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 321,1}, {-321,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 321,1}, {-321,1}};
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
	const Particle & pip = DS.decayProducts()[ix].at( sign*211)[0];
	const Particle & Kp  = DS.decayProducts()[ix].at( sign*321)[0];
	const Particle & Km  = DS.decayProducts()[ix].at(-sign*321)[0];
	double mplus  = (Kp.momentum()+pip.momentum()).mass2();
	double mminus = (Km.momentum()+pip.momentum()).mass2();
	double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	_h_KK[0]->fill(mKK);
	_h_KK[1]->fill(mKK);
	_h_Kppi->fill(mplus);
	_h_Kmpi->fill(mminus);
	_dalitz->fill(mKK,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK[0]);
      normalize(_h_KK[1],1.,false);
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK[2],_h_Kmpi,_h_Kppi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1830524);

}
