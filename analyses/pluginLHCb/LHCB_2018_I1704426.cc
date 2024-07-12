// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief 
  class LHCB_2018_I1704426 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2018_I1704426);


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
	book(_h[ix   ],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode =  { { 321,1}, { -321,1}, { 211,1}, { -211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = D0.decaying()[ix].pid()/421;
	if ( D0.modeMatches(ix,4,mode)) {
	  const Particles & Kp = D0.decayProducts()[ix].at( sign*321);
	  const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	  const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	  const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	  double mpipi = (pip[0].momentum()+pim[0].momentum()).mass();
	  if(mpipi>0.4802 && mpipi<0.5072) continue;
	  _h[0]->fill((Kp [0].momentum()+Km [0].momentum()).mass()/MeV);
	  _h[1]->fill(mpipi/MeV);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2018_I1704426);

}
