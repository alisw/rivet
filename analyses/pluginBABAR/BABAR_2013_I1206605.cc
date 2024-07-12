// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D+ -> K+K-pi+
  class BABAR_2013_I1206605 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2013_I1206605);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      DP.addStable(PID::K0S);
      DP.addStable(PID::ETA);
      DP.addStable(PID::ETAPRIME);
      declare(DP, "DP");
      // histos
      book(_h_Kppi,1,1,1);
      book(_h_KK  ,1,1,2);
      book(_h_Kmpi,1,1,3);
      book(_dalitz, "dalitz",50,0.,2.,50,0.9,3.1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-321,1}, { 211,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-321,1}, {-211,1}};
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      // loop over particles
      for(unsigned int ix=0;ix<DP.decaying().size();++ix) {
	int sign = 1;
	if (DP.decaying()[ix].pid()>0 && DP.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DP.decaying()[ix].pid()<0 && DP.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & Kp = DP.decayProducts()[ix].at( sign*321)[0];
	const Particle & Km = DP.decayProducts()[ix].at(-sign*321)[0];
	const Particle & pip= DP.decayProducts()[ix].at( sign*211)[0];
	double mminus = (Km.momentum()+pip.momentum() ).mass2();
	double mplus  = (Kp.momentum()+pip.momentum() ).mass2();
	double mKK    = (Kp.momentum()+Km.momentum()).mass2();
	_h_Kppi->fill(mplus);
	_h_Kmpi->fill(mminus);
	_h_KK  ->fill(mKK);
	_dalitz->fill(mminus,mKK); 
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppi);
      normalize(_h_Kmpi);
      normalize(_h_KK  );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppi,_h_Kmpi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2013_I1206605);

}
