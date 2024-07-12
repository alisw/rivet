// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi  -> pi+pi-pi0
  class KEDR_2022_I2593308 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KEDR_2022_I2593308);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "psi");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_pipi[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 111,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & pi0 = psi.decayProducts()[ix].at( 111);
	const Particles & pip = psi.decayProducts()[ix].at( 211);
	const Particles & pim = psi.decayProducts()[ix].at(-211);
	double mminus = (pim[0].momentum()+pi0[0].momentum()).mass();
	double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass();
	double mneut  = (pip[0].momentum()+pim[0].momentum()).mass();
	double Cminus = pim[0].momentum().p3().unit().dot(pi0[0].momentum().p3().unit());
	double Cplus  = pip[0].momentum().p3().unit().dot(pi0[0].momentum().p3().unit());
	double Cneut  = pip[0].momentum().p3().unit().dot(pim[0].momentum().p3().unit());
	if(Cminus+Cplus+Cneut>-1.075) continue; 
	if(Cneut >Cminus && Cneut >Cplus) _h_pipi[0]->fill(mneut );
	if(Cplus >Cminus && Cplus >Cneut) _h_pipi[1]->fill(mplus );
	if(Cminus>Cplus  && Cminus>Cneut) _h_pipi[2]->fill(mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h_pipi[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KEDR_2022_I2593308);

}
