// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K- pi+ pi0 (or CC)
  class BABAR_2006_I722905 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I722905);


    /// @name Analysis methods
    /// @{

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
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-321,1},{ 211,1}, {111,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-211,1}, {111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.modeMatches(ix,3,mode)) {
	  sign = 1;
	}
	else if  (D0.modeMatches(ix,3,modeCC)) {
	  sign =-1;
	}
	else
	  continue;
	const Particle & pi0= D0.decayProducts()[ix].at(111)[0];
	const Particle & Km = D0.decayProducts()[ix].at(-sign*321)[0];
	double m2 = (pi0.momentum()+Km.momentum()).mass2();
	if( (D0.decaying()[ix].pid()>0 && sign==  1) ||
	    (D0.decaying()[ix].pid()<0 && sign== -1) )
	  _h[0]->fill(m2);
	else
	  _h[1]->fill(m2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I722905);

}
