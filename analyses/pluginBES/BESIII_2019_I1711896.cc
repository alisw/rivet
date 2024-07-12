// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> Lambda eta pi+
  class BESIII_2019_I1711896 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1711896);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable(PID::ETA);
      LAMBDAC.addStable(PID::ETAPRIME);
      LAMBDAC.addStable(PID::LAMBDA);
      declare(LAMBDAC, "LAMBDAC");
      // histograms
      book(_h,1,1,1);
      book(_dalitz, "dalitz",50,1.5,3.5,50,2.5,5.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 221,1}, { PID::LAMBDA,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 221,1}, {-PID::LAMBDA,1}};
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for(unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
	int sign = 1;
	if (LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & lam = LAMBDAC.decayProducts()[ix].at( sign*PID::LAMBDA);
	const Particles & pip = LAMBDAC.decayProducts()[ix].at( sign*PID::PIPLUS);
	const Particles & eta = LAMBDAC.decayProducts()[ix].at(      PID::ETA);
	double mLamPi  = (lam[0].momentum()+pip[0].momentum()).mass2();
	double mLamEta = (lam[0].momentum()+eta[0].momentum()).mass2();
	_dalitz ->fill(mLamPi,mLamEta);
	_h->fill(sqrt(mLamPi));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1711896);

}
