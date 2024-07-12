#// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"
namespace Rivet {


  /// @brief D0 -> KS0 pi+ pi-
  class BELLE_2014_I1289224 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2014_I1289224);


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
      // Histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,0.3,3.2,50,0.3,3.2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // define the decay mode
      static const map<PdgId,unsigned int> & mode  = { { 310,1}, { 211,1},{-211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = D0.decaying()[ix].pid()/421;
	// KS0 pi+pi-
	if (!D0.modeMatches(ix,3,mode) ) continue;
	const Particle & pip= D0.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim= D0.decayProducts()[ix].at(-sign*211)[0];
	const Particle & K0 = D0.decayProducts()[ix].at( 310)[0];
	double mminus = (pim.momentum()+K0.momentum() ).mass2();
	double mplus  = (pip.momentum()+K0.momentum() ).mass2();
	double mpipi  = (pip.momentum()+pim.momentum()).mass2();
	_h[1]->fill(mplus);
	_h[2]->fill(mminus);
	_h[0]->fill(mpipi);
	_dalitz->fill(mminus,mplus); 
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2014_I1289224);

}
