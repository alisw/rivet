// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K-pi+pi+
  class FOCUS_2007_I750701 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2007_I750701);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid== 411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      declare(DP, "DP");
      // histos
      book(_h_low,1,1,1);
      book(_h_high,1,1,2);
      book(_dalitz, "dalitz",50,0.,3.1,50,0.3,3.1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2}, {-321,1} };
      static const map<PdgId,unsigned int> & modeCC = { {-211,2}, { 321,1} };
      // Loop over D+ mesons
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      for(unsigned int ix=0;ix<DP.decaying().size();++ix) {
	int sign = 1;
	if     ( DP.modeMatches(ix,3,mode  )) sign = 1;
	else if( DP.modeMatches(ix,3,modeCC)) sign =-1;
	else
	  continue;
	const Particles & pip = DP.decayProducts()[ix].at( sign*211);
	const Particle  & Km  = DP.decayProducts()[ix].at(-sign*321)[0];
	double mplus  = (Km.momentum() +pip[0].momentum()).mass2();
	double mminus = (Km.momentum() +pip[1].momentum()).mass2();
	if(mplus<mminus) swap(mplus,mminus);
	_h_low ->fill(mminus);
	_h_high->fill(mplus );
	_dalitz->fill(mminus,mplus );
	_dalitz->fill(mplus ,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_low );
      normalize(_h_high);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_high,_h_low;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2007_I750701);

}
