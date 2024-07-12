// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K_S0 K+ K-
  class BESIII_2020_I1799437 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1799437);


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
      book(_h_K0Km,1,1,3);
      book(_h_K0Kp,1,1,2);
      book(_h_KpKm,1,1,1);
      book(_dalitz, "dalitz",50,0.9,1.9,50,0.9,1.9);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-321,1}, { 310,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( ! D0.modeMatches(ix,3,mode)  ) continue;
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & K0 = D0.decayProducts()[ix].at(310);
	const Particles & Kp = D0.decayProducts()[ix].at( sign*321);
	const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	double mminus = (Km[0].momentum()+K0[0].momentum() ).mass2();
	double mplus  = (Kp[0].momentum()+K0[0].momentum() ).mass2();
	double mKK  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	_h_K0Kp->fill(mplus);
	_h_K0Km->fill(mminus);
	_h_KpKm->fill(mKK);
	_dalitz->fill(mKK,mplus); 
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_K0Km);
      normalize(_h_K0Kp);
      normalize(_h_KpKm);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0Km,_h_KpKm,_h_K0Kp;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1799437);

}
