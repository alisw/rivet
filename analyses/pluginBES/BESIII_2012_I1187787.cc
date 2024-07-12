// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta_c -> KS0 K+-pi-+
  class BESIII_2012_I1187787 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1187787);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==441);
      declare(ufs, "UFS");
      DecayedParticles etac(ufs);
      etac.addStable(PID::PI0);
      etac.addStable(PID::K0S);
      etac.addStable(PID::ETA);
      etac.addStable(PID::ETAPRIME);
      declare(etac, "etac");
      // histograms
      book(_h_K0pip,1,1,1);
      book(_h_K0Kp ,1,1,2);
      book(_h_Kppip,1,1,3);
      book(_dalitz, "dalitz",50,0.,8.,50,0.,8.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, {310,1}, {-211,1} };
      static const map<PdgId,unsigned int> & modeCC = { {-321,1}, {310,1}, { 211,1} };
      DecayedParticles etac = apply<DecayedParticles>(event, "etac");
      // loop over particles
      for(unsigned int ix=0;ix<etac.decaying().size();++ix) {
	int sign=1;
	if(etac.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if(etac.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & KS0 = etac.decayProducts()[ix].at(      310)[0];
	const Particle & pim = etac.decayProducts()[ix].at(-sign*211)[0];
	const Particle & Kp  = etac.decayProducts()[ix].at( sign*321)[0];
	double mplus  = (Kp .momentum() + pim.momentum()).mass2();
	double mminus = (KS0.momentum() + pim.momentum()).mass2();
	double mKK    = (Kp .momentum() + KS0.momentum()).mass2();
	_h_K0Kp ->fill(sqrt(mKK   ));
	_h_Kppip->fill(sqrt(mplus ));
	_h_K0pip->fill(sqrt(mminus));
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppip);
      normalize(_h_K0pip);
      normalize(_h_K0Kp );
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppip,_h_K0pip,_h_K0Kp;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1187787);

}
